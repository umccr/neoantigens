#!/usr/bin/env python

import sys
from os.path import dirname, abspath
import json
import csv
import requests
from ngs_utils import logger
import itertools
from pyensembl import EnsemblRelease
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

base = sys.argv[1]
pizzly_flat_filt_fpath = base + '-flat-filtered.tsv'
pizzly_json_fpath = base + '.json'
output_json = base + '.filt.json'
input_fasta = base + '.fusions.fasta'
output_fasta = base + '.fusions.filt.fasta'
output_bedpe_path = base + '-flat-filtered.bedpe'

ENSEMBL_VERSION = 75  # GRCh37. For GRCh38, please use 86
ebl = EnsemblRelease(ENSEMBL_VERSION)

# pizzly tsv:
"""
geneA.name  geneA.id         geneB.name  geneB.id         paircount  splitcount  transcripts.list
TFF1        ENSG00000160182  RPL7A       ENSG00000148303  2          4           ENST00000291527_0:551_ENST00000463740_29:1164;ENST00000291527_0:551_ENST00000323345_33:891
"""

# pizzly json:
"""
{
  "genes" : [
    {
      "geneA" : { "id" : "ENSG00000160182", "name" : "TFF1"},
      "geneB" : { "id" : "ENSG00000148303", "name" : "RPL7A"},
      "paircount" : 2,
      "splitcount" : 4,
      "transcripts" : [       # all combinations of transcripts
        {
          "fasta_record": "ENST00000291527_0:551_ENST00000463740_29:1164",
          "transcriptA": {"id" : "ENST00000291527", "startPos" : 0, "endPos" : 551, "edit" : -6, "strand" : true},
          "transcriptB": {"id" : "ENST00000463740", "startPos" : 29, "endPos" : 1164, "edit" : -4, "strand" : true},
          "support" : 6,
          "reads" : [0, 1, 2, 3, 4, 5]
        },
        {
          "fasta_record": "ENST00000291527_0:551_ENST00000323345_33:891",
          "transcriptA": {"id" : "ENST00000291527", "startPos" : 0, "endPos" : 551, "edit" : -6, "strand" : true},
          "transcriptB": {"id" : "ENST00000323345", "startPos" : 33, "endPos" : 891, "edit" : -4, "strand" : true},
          "support" : 6,
          "reads" : [0, 1, 2, 3, 4, 5]
        }
      ],
    },
    ...
  ]
}
"""

# output bedpe:
"""
#chr 5p  start 5p  end 5p    chr 3p  start 3p   end 3p    name of fusion    tier  strand 3p  strand 5p  quantitation
?        -1        ?         ?       -1         ?         TFF1>>RPL7A       -     ?          ?          -
"""

# Reading filtered tsv
filt_fusions = set()
with open(pizzly_flat_filt_fpath) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        filt_fusions.add((row['geneA.name'], row['geneB.name']))

""" Using pyensembl with a locally built database to map to genome coordinates
Please make sure to use the same Ensembl version as the one used with pizzly (75 for GRCh37, 86 for GRCh37). Same should go with INTEGRATE-Neo.
"""
class FusionPart:
    def __init__(self, transcript, t_start, t_end, g_start, g_end):
        self.transcript = transcript
        self.t_start = t_start
        self.t_end = t_end
        self.g_start = g_start
        self.g_end = g_end

    def fasta(self):
        return Seq(self.transcript.sequence[self.t_start:self.t_end])

def _tx_coord_to_genome_coord(transcript, coord):
    length = len(transcript)
    assert coord <= length, f'Coordinate {coord} > transcript length {length}, transcript: {transcript}'

    if transcript.strand == '-':
        coord = length - coord

    if not transcript.exons:
        logger.err(f'  No exons for transcript {transcript.transcript_id}')
        return None
    offset_remain = coord
    # print('  looking for coord', coord, f', in {len(transcript.exons)} exons, total length {length}')
    for exon in transcript.exons:
        # print('    exon len=', len(exon))
        # print('    offset_remain=', offset_remain)
        if offset_remain - len(exon) <= 0:
            # print('    returning exon.start + offset_remain = ', exon.start + offset_remain)
            return exon.start + offset_remain - 1   # -1 to convert from 1-based to 0-based
        offset_remain -= len(exon)
    assert False  # correct code should return

def _query_local_gtf(t_data):
    id = t_data['id']
    try:
        transcript = ebl.transcript_by_id(id)
    except:
        print(f'  Transcript {id} not found in Ensembl database')
        return None

    t_start = int(t_data['startPos'])
    t_end = int(t_data['endPos'])

    g_start = _tx_coord_to_genome_coord(transcript, t_start)
    g_end = _tx_coord_to_genome_coord(transcript, t_end)
    if not g_start:
        logger.err(f'  Error: could not find start coordinate {t_start} in transcript {id}')
        return None
    if not g_end:
        logger.err(f'  Error: could not find end coordinate {t_end} in transcript {id}')
        return None

    g_start, g_end = sorted([g_start, g_end])
    # print(f'  Transcript: {id} {transcript.strand}: {start}-{end} -> {g_chrom}:{g_start}-{g_end}')
    return FusionPart(transcript, t_start, t_end, g_start, g_end)

# Read json
filtered_data = []
with open(pizzly_json_fpath) as f:
    data = json.load(f)
    for g_event in data['genes']:
        gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
        if (gene_a, gene_b) in filt_fusions:
            filtered_data.append(g_event)

# Read fasta
fasta_dict = SeqIO.index(input_fasta, 'fasta')

# Write filtered json
with open(output_json, 'w') as f:
    json.dump({'genes': filtered_data}, f)

filt_fasta_records = []
split_fasta_records = []

# Write bedpe
with open(output_bedpe_path, 'w') as out:
    for g_event in filtered_data:
        gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
        print(gene_a + '>>' + gene_b)
        for t_event in g_event['transcripts']:
            t_a_data = t_event['transcriptA']  # {"id" : "ENST00000489283", "startPos" : 0, "endPos" : 168, "edit" : 0, "strand" : true}
            t_b_data = t_event['transcriptB']  # {"id" : "ENST00000333167", "startPos" : 496, "endPos" : 1785, "edit" : 7, "strand" : true}

            # a_chrom, a_start, a_end, a_strand = _query_coords(t_a_data)
            # b_chrom, b_start, b_end, b_strand = _query_coords(t_b_data)
            # print(f"  Transcript B from REST API: {t_b_data['id']}: {t_b_data['startPos']}-{t_b_data['endPos']} -> {b_chrom}:{b_start}-{b_end}, {b_strand}")
            a_part = _query_local_gtf(t_a_data)
            b_part = _query_local_gtf(t_b_data)
            if a_part and b_part:
                # print(f"  Transcript B from local DB: {t_b_data['id']}: {t_b_data['startPos']}-{t_b_data['endPos']} -> {b_chrom}:{b_start}-{b_end}, {b_strand}")

                bedpe_fields = a_part.transcript.contig, a_part.g_start, a_part.g_end, \
                               b_part.transcript.contig, b_part.g_start, b_part.g_end, \
                               gene_a + '>>' + gene_b, '-', a_part.transcript.strand, b_part.transcript.strand, '-'
                out.write('\t'.join(map(str, bedpe_fields)) + '\n')

                fasta_rec = fasta_dict[t_event['fasta_record']]
                filt_fasta_records.append(fasta_rec)

                # Split fasta records
                a_rec = SeqRecord(a_part.fasta(), f'{a_part.transcript.transcript_id}_{a_part.t_start}:{a_part.t_end}', '', '')
                split_fasta_records.append(a_rec)

                b_rec = SeqRecord(b_part.fasta(), f'{b_part.transcript.transcript_id}_{b_part.t_start}:{b_part.t_end}', '', '')
                split_fasta_records.append(b_rec)

                assert len(fasta_rec.seq) == len(a_rec.seq) + len(b_rec.seq)
                assert str(fasta_rec.seq) == str(a_rec.seq) + str(b_rec.seq), \
                    f'Seq {fasta_rec.id} != seq {a_rec.id} + seq {b_rec.id}: \n {fasta_rec.seq} \n != \n {a_rec.seq} \n + \n {b_rec.seq} \n full A transcript: \n {a_part.transcript.sequence} \n' \
                    f'Check that the Ensembl versions and genome builds match. You must use the same one as was run for pizzly.'

# Write split fasta
SeqIO.write(split_fasta_records, base + '.split.fasta', 'fasta')

# Write fasta
SeqIO.write(filt_fasta_records, output_fasta, 'fasta')








""" USING ENSEMBL REST API

We can use Ensembl REST API to map coordinates in a transcript to genomic coordinates in cDNA. 
Of course it means that coordinates spanning multiple transcripts will map to multiple segments in a genome, because it will need to jusp through introns.
For the second part of TFF1>>RPL7A fusion, pizzly reported 2 transcripts. They map to the following segments:

http://rest.ensembl.org/map/cdna/ENST00000463740/29..1164?content-type=application/json ->
{"mappings":[
    {"assembly_name":"GRCh38","end":133348246,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133348246},
    {"assembly_name":"GRCh38","end":133349042,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133348922},
    {"assembly_name":"GRCh38","end":133349700,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133349551},
    {"assembly_name":"GRCh38","end":133350052,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133349912},
    {"assembly_name":"GRCh38","end":133350727,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133350240},
    {"assembly_name":"GRCh38","end":133351071,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133351002},
    {"assembly_name":"GRCh38","end":133351426,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133351262}]}

http://rest.ensembl.org/map/cdna/ENST00000323345/33..891?content-type=application/json ->
{"mappings":[
    {"assembly_name":"GRCh38","end":133348246,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133348246},
    {"assembly_name":"GRCh38","end":133349042,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133348922},
    {"assembly_name":"GRCh38","end":133349700,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133349551},
    {"assembly_name":"GRCh38","end":133350052,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133349912},
    {"assembly_name":"GRCh38","end":133350319,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133350240},
    {"assembly_name":"GRCh38","end":133350727,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133350597},
    {"assembly_name":"GRCh38","end":133351071,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133351002},
    {"assembly_name":"GRCh38","end":133351426,"seq_region_name":"9","gap":0,"strand":1,"coord_system":"chromosome","rank":0,"start":133351262}]}

By playing around with the query, we see that we should take the coordinates from "start" of the first segment to the "end" of the last segment.
"""

url_tmpl = 'http://rest.ensembl.org/map/cdna/{tid}/{start}..{end}?'

def _query_coords(t_data):
    id = t_data['id']
    start = int(t_data['startPos']) + 1  # 0-based to 1-based
    end = int(t_data['endPos'])
    url = url_tmpl.format(tid=id, start=start, end=end)
    r = requests.get(url, headers={ "Content-Type" : "application/json"})

    g_chrom = None
    g_start = None
    g_end = None
    g_strand = None

    if not r.ok:
        r.raise_for_status()
        logger.err(f'Could not find transcript {id}')
    else:
        decoded = r.json()
        mappings = decoded['mappings']
        if not mappings:
            logger.err(f'No mappings for transcript {id} coordinates {start}:{end}')
        else:        
            g_chrom = mappings[0]['seq_region_name']
            g_start = int(mappings[0]['start']) - 1  # 1-based to 0-based
            g_end = int(mappings[-1]['end'])
            g_strand = '+' if mappings[0]['strand'] == '1' else '-'

    return g_chrom, g_start, g_end, g_strand
