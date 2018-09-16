#!/usr/bin/env python

import sys
from os.path import dirname, abspath, splitext
import json
import csv
import requests
from ngs_utils import logger
import itertools
from os import environ
from pyensembl import EnsemblRelease
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import math
import click


"""
Usage: pizzly_to_bedpe.py /path/to/bcbio/final/sample/pizzly/sample -o sample.bedpe

Assumes the following file to be under /path/to/bcbio/final/sample/pizzly/sample:

1. pizzly tsv `/path/to/bcbio/final/sample/pizzly/sample-flat-filtered.tsv`:

geneA.name  geneA.id         geneB.name  geneB.id         paircount  splitcount      transcripts.list
TFF1        ENSG00000160182  RPL7A       ENSG00000148303  2          4               ENST00000291527_0:551_ENST00000463740_29:1164;ENST00000291527_0:551_ENST00000323345_33:891

2. pizzly json `/path/to/bcbio/final/sample/pizzly/sample-flat-filtered.json`:

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

Generates on output `sample.bedpe`:

#chr 5p  start 5p  end 5p    chr 3p  start 3p   end 3p    name of fusion    tier  strand 3p  strand 5p  quantitation
21       -1        -1        9       -1         136215101 TFF1>>RPL7A	    -     -          +          -

"""

@click.command()
@click.argument('prefix')
@click.option('-o', '--output-bedpe', type=click.Path(), help='Output bedpe file path')
@click.option('--output-fasta', type=click.Path(), help='Filtered fasta file having only fusions from bedpe')
@click.option('--output-json', type=click.Path(), help='Filtered JSON file having only fusions from bedpe')
@click.option('-s', '--support', help='Minimal read support to keep an event', default=5)
@click.option('-g', '--genome', help='Genome build')
def main(prefix, output_bedpe, output_fasta=None, output_json=None, support=None, genome=None):
    pizzly_flat_filt_fpath = prefix + '-flat-filtered.tsv'
    pizzly_json_fpath = prefix + '.json'
    input_fasta = prefix + '.fusions.fasta'

    ebl = EnsemblRelease(75 if genome in ['GRCh37', 'hg19'] else 86)

    # Reading filtered tsv
    filt_fusions = set()
    with open(pizzly_flat_filt_fpath) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            filt_fusions.add((row['geneA.name'], row['geneB.name']))

    # Read json
    json_data = {'genes': []}
    with open(pizzly_json_fpath) as f:
        data = json.load(f)
        for g_event in data['genes']:
            gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
            if (gene_a, gene_b) in filt_fusions:
                json_data['genes'].append(g_event)

    # Read fasta
    fasta_dict = SeqIO.index(input_fasta, 'fasta')

    def _transcript_is_good(transcript):
        return transcript and transcript.biotype == 'protein_coding'

    filt_json_data = {'genes': []}
    filt_fasta_records = []
    filt_event_count = 0
    filt_transcript_event_count = 0

    # Write bedpe
    with open(output_bedpe, 'w') as bedpe_fh:
        for g_event in json_data['genes']:  # {'geneA', 'geneB', 'paircount', 'splitcount', 'transcripts', 'readpairs'}
            gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
            print(gene_a + '>>' + gene_b)

            # # first pass to select the longest transcripts
            # def _longest_tx(key):
            #     return max((ebl.transcript_by_id(te[f'transcript{key}']['id']) for te in g_event['transcripts']), key=lambda t: len(t))
            # a_tx = _longest_tx('A')
            # b_tx = _longest_tx('B')
            # print(f'Longest transcriptA: {a_tx.id}, Longest transcriptB: {b_tx.id}')
            # try:
            #     t_event = [te for te in g_event['transcripts'] if te['transcriptA']['id'] == a_tx.id and te['transcriptB']['id'] == b_tx.id][0]
            # except:
            #     print(f"No event with 2 longest transcripts. Available events: {', '.join(te['transcriptA']['id'] + 
            #           '>>' + te['transcriptB']['id'] for te in g_event['transcripts'])}")
            #     raise
                
            filt_g_event = {k: v for k, v in g_event.items() if k != 'readpairs'}
            filt_g_event['transcripts'] = []
            transcript_bedpe_entries = []  # collecting first to get rid of duplicates
            for t_event in g_event['transcripts']:
                t_a_data = t_event['transcriptA']  # {"id" : "ENST00000489283", "startPos" : 0, "endPos" : 168, "edit" : 0, "strand" : true}
                t_b_data = t_event['transcriptB']  # {"id" : "ENST00000333167", "startPos" : 496, "endPos" : 1785, "edit" : 7, "strand" : true}

                a_transcript = b_transcript = None
                try:
                    a_transcript = ebl.transcript_by_id(t_a_data['id'])
                except:
                    print(f"  Transcript A: {t_a_data['id']} not found in Ensembl database")
                try:
                    b_transcript = ebl.transcript_by_id(t_b_data['id'])
                except:
                    print(f"  Transcript B: {t_b_data['id']} not found in Ensembl database")

                if not _transcript_is_good(a_transcript) or not _transcript_is_good(b_transcript):
                    continue

                a_part = _query_local_gtf(a_transcript, int(t_a_data['startPos']), int(t_a_data['endPos']))
                b_part = _query_local_gtf(b_transcript, int(t_b_data['startPos']), int(t_b_data['endPos']))
                if not a_part or not b_part:
                    continue
                
                if t_event['support'] < support:
                    continue

                # for writing filtered json
                filt_g_event['transcripts'].append(t_event)
                filt_transcript_event_count += 1

                # writing bedpe
                bedpe_fields = [a_part.transcript.contig, a_part.g_start, a_part.g_end, \
                                b_part.transcript.contig, b_part.g_start, b_part.g_end, \
                                gene_a + '>>' + gene_b, 3, a_part.transcript.strand, b_part.transcript.strand, t_event['support']]
                # bedpe_fields += [a_part.transcript.transcript_id + ':' + str(len(a_part.transcript)), a_part.t_start, a_part.t_end]
                # bedpe_fields += [b_part.transcript.transcript_id + ':' + str(len(b_part.transcript)), b_part.t_start, b_part.t_end]
                bedpe_fields = tuple(bedpe_fields)
                if not bedpe_fields in transcript_bedpe_entries:  # ignoring duplicates
                    transcript_bedpe_entries.append(bedpe_fields)

                # for writing filtered fasta
                fasta_rec = fasta_dict[t_event['fasta_record']]
                filt_fasta_records.append(fasta_rec)

                # split fasta records to make sure we produce same the fasta with our coordinates as pizzly
                a_rec = SeqRecord(a_part.fasta(), f'{a_part.transcript.transcript_id}_{a_part.t_start}:{a_part.t_end}', '', '')
                b_rec = SeqRecord(b_part.fasta(), f'{b_part.transcript.transcript_id}_{b_part.t_start}:{b_part.t_end}', '', '')
                assert len(fasta_rec.seq) == len(a_rec.seq) + len(b_rec.seq)
                assert str(fasta_rec.seq) == str(a_rec.seq) + str(b_rec.seq), \
                    f'Seq {fasta_rec.id} != seq {a_rec.id} + seq {b_rec.id}: \n {fasta_rec.seq} \n != \n {a_rec.seq} \n \n ' \
                    f'{b_rec.seq} \n full A transcript: \n {a_part.transcript.sequence} \n' \
                    f'Check that the Ensembl versions and genome builds match. You must use the same one as was run for pizzly.'

            if not transcript_bedpe_entries:
                logger.warn(f'All transcript events filtered out for fusion {gene_a}>>{gene_b}, skipping')
            else:
                filt_json_data['genes'].append(filt_g_event)
                filt_event_count += 1
                for bedpe_fields in transcript_bedpe_entries:
                    bedpe_fh.write('\t'.join(map(str, bedpe_fields)) + '\n')

    # Write filtered json
    if output_json:
        with open(output_json, 'w') as f:
            json.dump(filt_json_data, f, indent=4)

    # Write fasta
    if output_fasta:
        SeqIO.write(filt_fasta_records, output_fasta, 'fasta')

    print()
    print(f'Written {filt_transcript_event_count} transcript events for {filt_event_count} fusions into bedpe: {output_bedpe}')


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
    if transcript.strand == '-':
        coord = length - coord

    if coord == 0:
        return -1
    if coord == length:
        return math.inf
    assert 0 < coord < length, f'Coordinate {coord} must be above 0 and below transcript length {length}, transcript: {transcript}'

    if not transcript.exons:
        logger.err(f'  No exons for transcript {transcript.transcript_id}')
        return None
    offset_remain = coord
    # print('  looking for coord', coord, f', in {len(transcript.exons)} exons, total length {length}')
    exons = transcript.exons
    if transcript.strand == '-':
        exons = reversed(exons)
    for exon in exons:
        assert offset_remain > 0
        # print('    exon len=', len(exon))
        # print('    offset_remain=', offset_remain)
        if offset_remain - len(exon) <= 0:
            # print('    returning exon.start + offset_remain = ', exon.start + offset_remain)
            return exon.start - 1 + offset_remain   # -1 to convert from 1-based to 0-based
        offset_remain -= len(exon)
    assert False  # correct code should return

def _query_local_gtf(ebl_transcript, t_start, t_end):
    # id = t_data['id']
    # try:
    # ebl_transcript = ebl_transcript or ebl.transcript_by_id(id)
    # except:
    #     print(f'  Transcript {id} not found in Ensembl database')
    #     return None

    # import pdb; pdb.set_trace()

    g_start = _tx_coord_to_genome_coord(ebl_transcript, t_start)
    g_end = _tx_coord_to_genome_coord(ebl_transcript, t_end)
    if g_start is None:
        logger.err(f'  Error: could not find start coordinate {t_start} in transcript {id}')
        return None
    if g_end is None:
        logger.err(f'  Error: could not find end coordinate {t_end} in transcript {id}')
        return None

    g_start, g_end = sorted([g_start, g_end])

    if g_end == math.inf:
        g_end = -1
    if g_start == g_end == -1:
        logger.warn(f'  Fusion takes the whole transcript {ebl_transcript.transcript_id} from the beginning to the end. ' +
                    f'That\'s suspicious, so we are skipping it.')
        return None

    # print(f'  Transcript: {id} {transcript.strand}: {start}-{end} -> {g_chrom}:{g_start}-{g_end}')
    return FusionPart(ebl_transcript, t_start, t_end, g_start, g_end)


if __name__ == '__main__':
    main()
