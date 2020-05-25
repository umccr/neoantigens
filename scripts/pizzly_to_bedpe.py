#!/usr/bin/env python
import os
import sys
from collections import defaultdict
from os.path import dirname, abspath, splitext, isfile, join
import json
import csv
import click
from ngs_utils import logger
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import safe_mkdir
from pyensembl import EnsemblRelease, Transcript
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
from memoized_property import memoized_property
from reference_data import api as refdata


"""
Map transcript coordinates to genomic coordinates, output a BEDPE file compatible with pVACfuse.

Please make sure to use the same Ensembl version as the one used with pizzly (75 for GRCh37, 86 for GRCh37). 
Same should go with INTEGRATE-Neo.

Usage: 
    pizzly_to_bedpe.py /path/to/bcbio/final/sample/pizzly/sample -o sample.bedpe

Assumes the following files to be under /path/to/bcbio/final/sample/pizzly/sample:

1. pizzly tsv 
```
cat /path/to/bcbio/final/sample/pizzly/sample-flat-filtered.tsv
geneA.name  geneA.id         geneB.name  geneB.id         paircount  splitcount      transcripts.list
TFF1        ENSG00000160182  RPL7A       ENSG00000148303  2          4               ENST00000291527_0:551_ENST00000463740_29:1164;ENST00000291527_0:551_ENST00000323345_33:891
```

2. pizzly json
```
cat /path/to/bcbio/final/sample/pizzly/sample-flat-filtered.json
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
```

Writes a bedpe file under the path provided with -o:

```
cat sample.bedpe
#chr 5p  start 5p  end 5p    chr 3p  start 3p   end 3p    name of fusion    tier  strand 3p  strand 5p  quantitation
21       -1        -1        9       -1         136215101 TFF1>>RPL7A	    -     -          +          -
```
"""

ENSEMBL_RELEASE=95

@click.command()
@click.argument('prefix')
@click.option('-o', '--output-bedpe', type=click.Path(), help='Output bedpe file path')
@click.option('--output-fasta', type=click.Path(), help='Filtered fasta file having only fusions from bedpe')
@click.option('--output-json', type=click.Path(), help='Filtered JSON file having only fusions from bedpe')
@click.option('-s', '--min-read-support', help='Minimal read support to keep an event', default=5)
@click.option('-e', '--ensembl-release', default=91, help='Set to 75 for GRCh37')
@click.option('-p', '--peptide-flanking-len', type=int,
              help='Reported only a part of fusion peptide around the breakpoint. Take `p` number of aminoacids '
                   'from each side of the fusion. Plus the junction peptide (the resulting peptide length will be `p`*2+1).')
@click.option('-d', '--debug', is_flag=True)
@click.option('--no-filtering', 'no_filtering', is_flag=True)
@click.option('--transcript-check/--no-transcript-check', 'check_transcript', is_flag=True, default=True)

@click.option('--trx-fa', type=click.Path(), help='Reference transcriptome fasta')
@click.option('-r', '--reads', type=click.Path(), help='Fastq files with reads for re-quantifying', multiple=True)
@click.option('-t', '--min-tpm', help='Minimal TPM for a transcript. Re-quanifies with kallisto, '
                                      'thus applies only if reads are available.', default=1)

# @click.option('--keep-noncoding', is_flag=True,
#               help='Keep fusions that do not produce a peptide around the junction')
def main(prefix, output_bedpe, output_fasta=None, output_json=None, min_read_support=None, ensembl_release=None,
         peptide_flanking_len=None, debug=False, no_filtering=False, check_transcript=True,
         trx_fa=None, reads=None, min_tpm=None):

    # input_flat_fpath = prefix + '-flat.tsv'
    input_json_fpath = prefix + '.json'
    input_fasta = prefix + '.fusions.fasta'
    output_bedpe = abspath(output_bedpe)

    logger.init(debug)

    global ENSEMBL_RELEASE
    ENSEMBL_RELEASE = ensembl_release
    ebl = EnsemblRelease(ENSEMBL_RELEASE)

    # Reading filtered tsv
    # filt_fusions = set()
    # with open(input_flat_fpath) as f:
    #     for row in csv.DictReader(f, delimiter='\t'):
    #         filt_fusions.add((row['geneA.name'], row['geneB.name']))

    # Read json
    json_data = {'genes': []}
    with open(input_json_fpath) as f:
        data = json.load(f)
        for g_event in data['genes']:
            gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
            # if (gene_a, gene_b) in filt_fusions:
            json_data['genes'].append(g_event)

    # Read fasta
    fasta_dict = SeqIO.index(input_fasta, 'fasta')

    # First round: genomic coordinates and fasta
    logger.info(f'Round 1: reading {len(json_data["genes"])} gene-pairs events from pizzly JSON')
    fusions = []
    for g_event in json_data['genes']:  # {'geneA', 'geneB', 'paircount', 'splitcount', 'transcripts', 'readpairs'}
        gene_a, gene_b = g_event['geneA']['name'], g_event['geneB']['name']
        # logger.info(f'Processing event {gene_a}>>{gene_b}')

        met_fasta_keys = set()  # collecting to get rid of duplicate transcript events
        for t_event in g_event['transcripts']:
            fusion = Fusion.create_from_pizzly_event(ebl, t_event)

            if check_transcript:
                if not _transcript_is_good(fusion.side_5p.trx) or not _transcript_is_good(fusion.side_3p.trx):
                    logger.info(f'Transcripts {fusion.side_5p.trx} and {fusion.side_3p.trx} didn\'t pass check')
                    continue

            if no_filtering is not True and fusion.support < min_read_support: continue

            if not fusion.calc_genomic_positions(): continue

            # comparing our fasta to pizzly fasta
            fusion.fasta_rec = fasta_dict[t_event['fasta_record']]
            _check_fusion_fasta(fusion.fasta_rec, fusion)

            # skipping duplicate fastas
            k = fusion.side_5p.trx.id, fusion.side_3p.trx.id, fusion.fasta
            assert k not in met_fasta_keys
            met_fasta_keys.add(k)

            fusions.append(fusion)
        # if not met_fasta_keys:
        #     logger.info('   Filtered all fusions for this gene pair.')
        if met_fasta_keys:
            logger.info(f'Keeping {len(met_fasta_keys)} fusion(s)')

    if not fusions:
        logger.warn('Finished: no fusions passed filtering')
        sys.exit(0)

    # Calculate expression of fused transcripts
    expr_by_fusion = None
    if reads:
        # filtered fasta for re-calling expression
        work_dir = safe_mkdir(splitext(output_bedpe)[0] + '_quant')
        fasta_path = join(work_dir, 'fusions.fasta')
        fasta_recs = [f.fasta_rec for f in fusions]
        SeqIO.write(fasta_recs, fasta_path, 'fasta')

        trx_fa = trx_fa or splitext(get_ref_file('hg38', 'gtf'))[0] + '.fa'
        assert isfile(trx_fa)
        expr_by_fusion = requanitify_pizzly(trx_fa, fasta_path, work_dir, reads)
        # expr_by_fusion = {fusion-fasta-id -> {length  eff_length  est_counts   tpm}}

    # Second round: peptides and expression
    logger.info()
    logger.info(f'Round 2: making peptides for {len(fusions)} events in '
                f'{len(set([(f.side_3p.trx.gene.name, f.side_5p.trx.gene.name) for f in fusions]))} genes pairs')
    met_peptide_keys = set()  # collecting to get rid of duplicate peptides
    bedpe_entries = []
    peptide_fusions = []
    for fusion in fusions:
        if fusion.side_3p.trx.contains_start_codon:
            logger.info(f'Translating {fusion.side_5p.trx.gene.name}>>{fusion.side_3p.trx.gene.name} fusion: {fusion}')
            fusion.make_peptide(peptide_flanking_len)
            if fusion.peptide:
                _verify_peptides(fusion.fasta_rec, fusion, peptide_flanking_len)

            # skipping duplicate peptides
            k = fusion.side_5p.trx.gene.name, fusion.side_3p.trx.gene.name, fusion.peptide
            if k in met_peptide_keys: continue
            met_peptide_keys.add(k)

        # writing bedpe
        entry = fusion.to_bedpe()

        # add expression
        if expr_by_fusion:
            entry.update(expr_by_fusion[fusion.fasta_rec.id])
            if no_filtering is not True and float(entry['tpm']) < min_tpm: continue

        if fusion.peptide:
            peptide_fusions.append(fusion)
        bedpe_entries.append(entry)

    # Writing bedpe
    with open(output_bedpe, 'w') as bedpe_fh:
        bedpe_header = [
            'chr 5p',
            'start 5p',
            'end 5p',
            'chr 3p',
            'start 3p',
            'end 3p',
            'name',
            'tier',
            'strand 5p',
            'strand 3p',
            'support',
            'is canon bndry',
            'inframe',
            'peptide',
            'fusion pos',
            'nt in the break',
            'transcripts',
            'is canon intron dinuc',
        ]
        if expr_by_fusion:
            bedpe_header.extend(list(expr_by_fusion.values())[0].keys())
        bedpe_writer = csv.DictWriter(bedpe_fh, fieldnames=bedpe_header, delimiter='\t')
        bedpe_writer.writeheader()
        for bedpe_entry in bedpe_entries:
            bedpe_writer.writerow(bedpe_entry)

    # _test_pvac(output_bedpe)

    # Write fasta
    if output_fasta:
        SeqIO.write([f.fasta_rec for f in peptide_fusions], output_fasta, 'fasta')

    logger.info()
    logger.info(f'Written {len(peptide_fusions)} fusions in '
                f'{len(set([(f.side_3p.trx.gene.name, f.side_5p.trx.gene.name) for f in peptide_fusions]))} '
                f'gene pairs good peptides bedpe: {output_bedpe}')


def _test_pvac(bedpe_path):
    pvac_bedpe = bedpe_path.replace('.bedpe', '.pvac.bedpe')
    pvac_tsv_path = bedpe_path.replace('.bedpe', '.pvac.tsv')
    pvac_fasta_fpath = bedpe_path.replace('.bedpe', '.pvac.fasta')
    pvac_fasta_key_fpath = bedpe_path.replace('.bedpe', '.pvac.fasta_key')

    run_simple(f'grep -v ^chr {bedpe_path} > {pvac_bedpe}')

    from lib.fasta_generator import FusionFastaGenerator
    from lib.pipeline import MHCIPipeline

    class_i_arguments = {
        'input_file':             pvac_bedpe,
        'input_file_type':       'bedpe',
        'sample_name':            bedpe_path.replace('.bedpe', '.pvac'),
        'alleles':               'HLA-A*02:01',
        'prediction_algorithms': 'NetMHCcons',
        'output_dir':             dirname(bedpe_path),
        'epitope_lengths':        11,
    }
    if isfile(pvac_tsv_path): os.remove(pvac_tsv_path)
    pipeline = MHCIPipeline(**class_i_arguments)
    pipeline.convert_vcf()

    generate_fasta_params = {
        'input_file'                : pvac_tsv_path,
        'epitope_length'            : 11,
        'output_file'               : pvac_fasta_fpath,
        'output_key_file'           : pvac_fasta_key_fpath,
        'downstream_sequence_length': 1000,
    }
    fasta_generator = FusionFastaGenerator(**generate_fasta_params)
    fasta_generator.execute()

    generate_fasta_params = {
        'input_file'                : pvac_tsv_path,
        'epitope_length'            : 8,
        'output_file'               : pvac_fasta_fpath,
        'output_key_file'           : pvac_fasta_key_fpath,
        'downstream_sequence_length': 1000,
    }
    fasta_generator = FusionFastaGenerator(**generate_fasta_params)
    fasta_generator.execute()


def _transcript_is_good(trx):
    return \
        trx is not None and \
        (trx.biotype == 'protein_coding' or 'decay' in trx.biotype) and \
        trx.support_level == 1


def _find_transcript(ebl, id, name):
    return ebl.transcript_by_id(id)
    # try:
    #     return ebl.transcript_by_id(id)
    # except:
    #     traceback.print_exc()
    #     logger.critical(f'Warning: transcript {name}: {id} not found in Ensembl database')


class FusionSide:
    """ Represents a one side of a fusion (belinging to one of 2 transcripts being fused)
    """
    def __init__(self, transcript: Transcript, bp_offset: int):
        self.trx = transcript
        self.bp_offset = bp_offset

        self.bp_genomic_pos = None
        self.bp_is_in_intron = None

    def calc_genomic_bp_pos(self):
        genomic_coord, is_in_intron = FusionSide.offset_to_genome_coord(self.trx, self.bp_offset)
        if genomic_coord is None:
            logger.critical(f'  Error: could not convert transcript {id} offest {genomic_coord} to genomic coordinate')
            return False

        if genomic_coord == -1:
            logger.warn(f'  Fusion in {self} takes the whole transcript {self.trx.id}. That\'s suspicious, so we are skipping it.')
            return False

        self.bp_genomic_pos = genomic_coord
        self.bp_is_in_intron = is_in_intron
        return True

    @staticmethod
    def offset_to_genome_coord(trx, offset):
        genomic_coord = None
        is_in_intron = None

        length = len(trx)
        offset = offset if trx.strand == '+' else length - offset
        if offset == 0 or offset == length:
            return -1, None

        assert 0 < offset < length, f'Coordinate {offset} must be above 0 and below transcript length {length}, ' \
                                    f'transcript: {trx}'
        if not trx.exons:
            logger.err(f'  No exons for transcript {trx.id}')
            return None, None

        offset_remain = offset
        # print('  looking for coord', coord, f', in {len(transcript.exons)} exons, total length {length}')
        exons = trx.exons
        if trx.strand == '-':
            exons = reversed(exons)
        for exon in exons:
            assert offset_remain > 0
            # print('    exon len=', len(exon))
            # print('    offset_remain=', offset_remain)
            next_offset_remain = offset_remain - len(exon)
            if next_offset_remain <= 0:
                # print('    returning exon.start + offset_remain = ', exon.start + offset_remain)
                genomic_coord = exon.start - 1 + offset_remain   # -1 to convert from 1-based to 0-based
                is_in_intron = next_offset_remain == 0
                break
            offset_remain = next_offset_remain
        assert genomic_coord is not None  # correct code should always produce something
        return genomic_coord, is_in_intron

    @staticmethod
    def create_5p(ebl_db, t_event):
        t_data = t_event['transcriptA']  # {"id" : "ENST00000489283", "startPos" : 0, "endPos" : 168, "edit" : 0, "strand" : true}
        trx = ebl_db.transcript_by_id(t_data['id'])

        start = int(t_data['startPos'])
        end = int(t_data['endPos'])

        assert start == 0, f'For 5\' transcript {trx.id}, start pos = {start} != 0. '
        return FusionSide(trx, bp_offset=end)

    @staticmethod
    def create_3p(ebl_db, t_event):
        t_data = t_event['transcriptB']  # {"id" : "ENST00000333167", "startPos" : 496, "endPos" : 1785, "edit" : 7, "strand" : true}
        trx = ebl_db.transcript_by_id(t_data['id'])

        start = int(t_data['startPos'])
        end = int(t_data['endPos'])

        assert end == len(trx), f'The transcript {trx.id} length in pizzly results ({end}) mismatches this transcript ' \
                                f'length {len(trx)} in Ensembl release {ENSEMBL_RELEASE}. ' \
                                f'Double check if the same Ensembl release is used in pizzly.'
        return FusionSide(trx, bp_offset=start)

class Fusion:
    def __init__(self, side_5p: FusionSide, side_3p: FusionSide, support: int, tier=3):
        self.side_5p = side_5p  # 5' transcript
        self.side_3p = side_3p  # 3' transcript

        self.support = support  # number of reads supporting the variant
        self.tier = tier

        self.is_canonical_boundary = None
        self.is_inframe = None
        self.peptide = None
        self.fusion_offset_in_peptide = None
        self.num_of_nt_in_the_break = None

        self.fasta_rec = None

    @staticmethod
    def create_from_pizzly_event(ebl_db, t_event):
        side_5p = FusionSide.create_5p(ebl_db, t_event)
        side_3p = FusionSide.create_3p(ebl_db, t_event)
        return Fusion(side_5p, side_3p, t_event['support'])

    def calc_genomic_positions(self):
        return self.side_3p.calc_genomic_bp_pos() and self.side_5p.calc_genomic_bp_pos()

    def to_bedpe(self):
        self.is_canonical_boundary = self.side_5p.bp_is_in_intron and self.side_3p.bp_is_in_intron
        assert self.side_5p.bp_genomic_pos
        assert self.side_3p.bp_genomic_pos

        entry = {
            'chr 5p':                self.side_5p.trx.contig,
            'start 5p':              -1 if self.side_5p.trx.strand == '+' else self.side_5p.bp_genomic_pos,
            'end 5p':                -1 if self.side_5p.trx.strand == '-' else self.side_5p.bp_genomic_pos,
            'chr 3p':                self.side_3p.trx.contig,
            'start 3p':              self.side_3p.bp_genomic_pos if self.side_3p.trx.strand == '+' else -1,
            'end 3p':                self.side_3p.bp_genomic_pos if self.side_3p.trx.strand == '-' else -1,
            'name':                  self.side_5p.trx.gene.name + '>>' + self.side_3p.trx.gene.name,
            'tier':                  self.tier,
            'strand 5p':             self.side_5p.trx.strand,
            'strand 3p':             self.side_3p.trx.strand,
            'support':               self.support,
            'is canon bndry':        'NA',
            'inframe':               'NA',
            'peptide':               'NA',
            'fusion pos':            'NA',
            'nt in the break':       'NA',
            'transcripts':           'NA',
            'is canon intron dinuc': 'NA',
        }
        if self.peptide:
            # ENST00000304636|ENST00000317840;ENST00000377795;ENST00000009530|ENST00000353334
            # 5' transcripts                 ;3' transcripts ;3' frameshift transcripts
            trx_line = self.side_5p.trx.transcript_id + ';' + \
                       (self.side_3p.trx.id if self.is_inframe else '') + ';' + \
                       (self.side_3p.trx.id if not self.is_inframe else '')

            entry.update({
                'is canon bndry':        '1' if self.is_canonical_boundary else '0',
                'inframe':               '1' if self.is_inframe else '0',
                'peptide':               self.peptide,
                'fusion pos':            self.fusion_offset_in_peptide,
                'nt in the break':       self.num_of_nt_in_the_break,
                'transcripts':           trx_line,
            })
        # fields += [self.side_a.transcript.transcript_id + ':' + str(len(self.side_a.transcript)),
        #            self.side_a.t_start, self.side_a.t_end]
        # fields += [self.side_b.transcript.transcript_id + ':' + str(len(self.side_b.transcript)),
        #            self.side_b.t_start, self.side_b.t_end]
        return entry

    @memoized_property
    def fasta(self):
        seq = self.side_5p.trx.sequence[:self.side_5p.bp_offset] + self.side_3p.trx.sequence[self.side_3p.bp_offset:]
        return Seq(seq, IUPAC.unambiguous_dna)

    def make_peptide(self, peptide_flanking_len=None):
        # 5' fasta
        if not self.side_5p.trx.contains_start_codon:
            logger.warn('No start codong in 5\' transcript')
            return

        transl_start = self.side_5p.trx.first_start_codon_spliced_offset
        if self.side_5p.bp_offset < transl_start:  # if the bp (t_end) falls before the beginning of translation
            return
        cds_seq_5p = Seq(self.side_5p.trx.sequence[transl_start:self.side_5p.bp_offset])
        fs_5p = len(cds_seq_5p) % 3
        if fs_5p != 0: logger.debug(f'  Frameshift of 5p sequence: {fs_5p}')

        pep_5p = _translate_from_start_codon(cds_seq_5p, to_stop=False, name='5\' fasta')
        if '*' in pep_5p:
            logger.info(f'   5\' petide has a STOP codon before breakpoint. Skipping.')
            assert min(self.side_5p.trx.stop_codon_spliced_offsets) < self.side_5p.bp_offset, \
                'We also expect pyensembl to report a STOP codon before the breakpoint'
            return

        # 3' fasta. Getting full sequence in case if it's an FS event that will produce a novel stop codon
        fs_3p = (self.side_3p.bp_offset - self.side_3p.trx.first_start_codon_spliced_offset) % 3
        if fs_3p != 0: logger.debug(f'  Frameshift of 3p sequence: {fs_3p}')
        seq_3p = Seq(self.side_3p.trx.sequence[self.side_3p.bp_offset:])

        # checking if the fusion produced a frameshift
        fusion_fs = (fs_5p + fs_3p) % 3
        if fusion_fs != 0: logger.debug(f"  Result fusion frameshift:  {fusion_fs}")
        is_inframe = fusion_fs == 0

        # junction peptide
        junction_codon = cds_seq_5p[len(cds_seq_5p)-fs_5p:]
        if junction_codon:  # == fs_5p != 0:
            start_3p_from = 3 - fs_5p
            junction_codon += seq_3p[:start_3p_from]
            junction_pep = junction_codon.translate()
            if junction_pep == '*':
                logger.info(f'   Junction codon is STOP, skipping')
                return
        else:
            junction_pep = ''
            start_3p_from = 0

        # 3' peptide
        pep_3p = _trim3(seq_3p[start_3p_from:]).translate()
        if pep_3p[0] == '*':
            logger.info(f'   The new 3\' peptide starts from STOP ({seq_3p[start_3p_from:][:3]} '
                        f'at position {self.side_3p.bp_offset}+{start_3p_from}), skipping translation.')
            return
        if '*' not in pep_3p:
            logger.info(f'   No STOP codon in fused peptide, skipping translation.')
            return
        pep_3p = _trim3(seq_3p[start_3p_from:]).translate(to_stop=True)

        logger.debug(f'  5\' peptide (len={len(pep_5p)}): '
                     f'{pep_5p if len(pep_5p) < 99 else pep_5p[:48] + "..." + pep_5p[-48:]}')
        if junction_pep:
            logger.debug(f'  Junction peptide: {junction_pep}')
        logger.debug(f'  3\' peptide{f" (shifted by {fusion_fs} from original)" if fusion_fs else ""} (len={len(pep_3p)}): '
                     f'{pep_3p if len(pep_3p) < 99 else pep_3p[:48] + "..." + pep_3p[-48:]}')

        # fusion peptide
        if peptide_flanking_len:
            # taking $(peptide_chunk_len) aminoacids from 5':
            pep_5p = pep_5p[-peptide_flanking_len:]
            # trying to make the total peptide to be of length $(peptide_chunk_len)*2+1:
            pep_3p = pep_3p[:peptide_flanking_len + 1 - len(junction_pep)] if is_inframe else pep_3p

        fusion_pep = pep_5p + junction_pep + pep_3p
        assert '*' not in fusion_pep

        self.peptide = fusion_pep
        self.is_inframe = is_inframe
        self.fusion_offset_in_peptide = peptide_flanking_len or len(pep_5p)
        self.num_of_nt_in_the_break = fs_5p

    def __repr__(self):
        return f'Fusion(' \
               f'{self.side_5p.trx.id}({self.side_5p.trx.strand}):{self.side_5p.bp_offset} >> ' \
               f'{self.side_3p.trx.id}({self.side_3p.trx.strand}):{self.side_3p.bp_offset})'


def _check_fusion_fasta(pizzly_fasta_rec, fusion):
    # split fasta records to make sure we produce same the fasta with our coordinates as pizzly
    fus_fasta_rec = SeqRecord(
        fusion.fasta,
        f'{fusion.side_5p.trx.id}_{fusion.side_5p.bp_offset}>>{fusion.side_3p.trx.id}_{fusion.side_3p.bp_offset}',
        '', '')
    # a_rec = SeqRecord(fusion.side_a.fasta,
    #                   f'{fusion.side_a.transcript.transcript_id}_{fusion.side_a.t_start}:{fusion.side_a.t_end}',
    #                   '', '')
    # b_rec = SeqRecord(fusion.side_b.fasta,
    #                   f'{fusion.side_b.transcript.transcript_id}_{fusion.side_b.t_start}:{fusion.side_b.t_end}',
    #                   '', '')
    assert len(pizzly_fasta_rec.seq) == len(fus_fasta_rec.seq)
    assert str(pizzly_fasta_rec.seq) == str(fus_fasta_rec.seq), \
        f'Seq {pizzly_fasta_rec.id} != seq {fusion.side_5p.trx.id} + seq {fusion.side_3p.trx.id}: \n ' \
        f'{pizzly_fasta_rec.seq} \n != \n {fus_fasta_rec.seq} \n ' \
        f'Check that the Ensembl versions and genome builds are matching. ' \
        f'You must use the same one as was run for pizzly.'


# trim seq to be a length of a multiple of 3 - to aboid BioPython translation warning about trailing nt
def _trim3(seq):
    return seq[:len(seq) // 3 * 3]


def _translate_from_start_codon(seq, to_stop, name):
    """ Seq must start with START. Translates until STOP.
    """
    codon_table = CodonTable.unambiguous_dna_by_name['Standard']
    if str(seq[:3]).upper() not in codon_table.start_codons:
        logger.critical(name + ' expected to start with a START codon: ' + seq[:3])
    pep_5p = _trim3(seq).translate(to_stop=to_stop)
    # for the case if the peptide starts with an alternative start codon, replace it with M
    return 'M' + pep_5p[1:]


def _verify_peptides(pizzly_fasta_rec, fusion, peptide_flanking_len=None):
    """ We can also calculate peptides directly from Pizzly fasta. Making sure we would produce identical peptides.
        In the future, we can ommit other logic and use Pizzly fasta directly.
    """
    pizzly_seq = pizzly_fasta_rec.seq
    start_codon_offset = fusion.side_5p.trx.first_start_codon_spliced_offset
    logger.debug(f'start_codon_offset: {start_codon_offset}')
    pizzly_seq = pizzly_seq[start_codon_offset:]
    bp_offset = fusion.side_5p.bp_offset - start_codon_offset
    logger.debug(f'bp_offset from start codon: {bp_offset}')
    if bp_offset <= 0:
        return

    pep = _translate_from_start_codon(pizzly_seq, to_stop=True, name='Pizzly fasta')
    assert pep[0] == 'M'

    fus_pos = bp_offset // 3
    if len(pep) < fus_pos:
        return

    fs_5p = bp_offset % 3
    fs_3p = (fusion.side_3p.bp_offset - fusion.side_3p.trx.first_start_codon_spliced_offset) % 3
    fusion_fs = (fs_5p + fs_3p) % 3
    if fusion_fs != 0: logger.debug(f"  Result fusion frameshift:  {fusion_fs}")
    is_inframe = fusion_fs == 0
    assert is_inframe == fusion.is_inframe, (is_inframe, fusion.is_inframe)

    if peptide_flanking_len:
        if is_inframe:
            logger.debug(f'{fusion} is in frame')
            pep = pep[fus_pos-peptide_flanking_len:fus_pos+peptide_flanking_len+1]
        else:
            logger.debug(f'{fusion} is off-frame')
            pep = pep[fus_pos-peptide_flanking_len:]
        fus_pos = peptide_flanking_len
        logger.debug(f'fusion position in trimmed peptide: {fus_pos}')

    assert fus_pos == fusion.fusion_offset_in_peptide, (fus_pos, fusion.fusion_offset_in_peptide)
    assert fusion.peptide == pep, (fusion.peptide, pep, fusion)


def requanitify_pizzly(ref_fa, fusions_fasta, work_dir, fastq):
    """ Returns dict fusion-fasta-id -> {length  eff_length  est_counts   tpm}
    """
    trx_with_fusions = join(work_dir, 'transcripts_with_fusions.fasta.gz')
    kidx = join(work_dir, 'transcripts_with_fusions.kidx')

    if not isfile(trx_with_fusions):
        run_simple(f"cat {ref_fa} {fusions_fasta} | gzip -c > {trx_with_fusions}")

    if not isfile(kidx):
        run_simple(f"kallisto index -k31 -i {kidx} {trx_with_fusions}")

    abundance = join(work_dir, 'abundance.tsv')
    if not isfile(abundance):
        run_simple(f"kallisto quant -i {kidx} -o {work_dir} {' '.join(fastq)}")

    logger.debug(f'Reading expression from {abundance}')
    expr_by_fusion = dict()
    with open(abundance) as f:
        header = f.readline().strip().split('\t')
        for row in csv.DictReader(f, delimiter='\t', fieldnames=header):
            expr_by_fusion[row['target_id']] = row
    return expr_by_fusion


if __name__ == '__main__':
    main()







