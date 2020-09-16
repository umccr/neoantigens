from os.path import isfile, join, dirname, abspath
import csv


def package_path():
    return dirname(abspath(__file__))


def nag_summarize(sample_names, quality_fpaths, output_fpath):
    cols = ['Sample', 'NeoantigenID', 'MT.Peptide.Form', 'NeoantigenQuality',
            'NeoantigenAlignment', 'IEDB_EpitopeAlignment', 'AlignmentScore', 'IEDB_Epitope']

    top_epitope_by_sample = dict()
    for sn, inp_fpath in zip(sample_names, quality_fpaths):
        with open(inp_fpath) as f:
            csv_reader = csv.DictReader(f, delimiter='\t')
            epitopes = []
            for e in csv_reader:
                epitopes.append(e)
            if epitopes:
                epitopes.sort(key=lambda _e: _e['NeoantigenQuality'])
                best = epitopes[-1]  # the top one?
            else:
                best = {c: '.' for c in cols}
                best['Sample'] = sn
            top_epitope_by_sample[best['Sample']] = best

    assert len(top_epitope_by_sample) == len(sample_names)
    assert len(top_epitope_by_sample) > 0, top_epitope_by_sample
    with open(output_fpath, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for sn in sample_names:
            out.write('\t'.join(top_epitope_by_sample[sn][c] for c in cols) + '\n')


