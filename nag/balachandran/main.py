#!/usr/bin/env python
import csv
import math
from os.path import join
import click
import pandas as pd
from Aligner import Aligner


@click.command()
@click.option('-l', '--list', 'epitope_list', help='Epitope list (columns: )')
@click.option('-a', '--alignments', 'alignments_dir', help='Folder with precomputed blast alignments')
@click.option('-A', '--alignment-score-threshold', 'alignment_score_threshold', type=click.FLOAT,
              default=26.0, help='Midpoint parameter of the logistic function, alignment score threshold')
@click.option('-K', '--slope-parameter', 'slope_parameter', type=click.FLOAT,
              default=1.0, help='Slope parameter of the logistic function')
@click.option('-o', '--output-file', 'output_file', help='Output file')

def main(epitope_list=None, alignments_dir=None, alignment_score_threshold=None, slope_parameter=None,
         output_file=None):

    # Compute MHC amplitudes for all neoantigens
    a_val_by_index = {}
    peptide_by_index = {}
    sample_by_index = {}

    with open(epitope_list) as f:
        for data in csv.DictReader(f, delimiter='\t'):
            index = data['id']
            sample = data['sample']
            mtpeptide = data['epitope']
            kdwt = data['wt_score']
            kdmt = data['mt_score']
            kdmt = float(kdmt)
            if kdwt == 'nan':
                kdwt = 1000.
            kdwt = float(kdwt)
            index = int(index)
            peptide_by_index[index] = mtpeptide.upper()
            a_val_by_index[index] = kdwt / kdmt
            sample_by_index[index] = sample

    # Compute TCR-recognition probabilities for all neoantigens
    aligner = Aligner()
    for sname in set(sample_by_index.values()):
        xml_path = join(alignments_dir, f'neoantigens_{sname}_iedb.xml')
        aligner.read_all_blast_alignments(xml_path)
    aligner.compute_rval(alignment_score_threshold, slope_parameter)

    # Compute qualities for all epitopes and write the result
    with open(output_file, 'w') as out:
        header = ['Sample', 'NeoantigenID', 'MT.Peptide.Form', 'NeoantigenQuality',
                  'NeoantigenAlignment', 'IEDB_EpitopeAlignment', 'AlignmentScore', 'IEDB_Epitope']
        out.write('\t'.join(header) + '\n')
        for index, peptide in peptide_by_index.items():
            a_val = a_val_by_index[index]
            [r_val, species, alignment] = aligner.get_rval(index)

            neo_alignment = alignment[0]
            epitope_alignment = alignment[1]
            score = alignment[2]

            quality = a_val * r_val
            res = [sample_by_index[index], index, peptide, quality, neo_alignment, epitope_alignment, score, species]
            out.write('\t'.join(map(str, res)) + '\n')


if __name__ == '__main__':
    main()
