"""
##################
## Installation ##
conda create -n pvac -c vladsaveliev -с bioconda -с conda-forge \
    snakemake-minimal vcfstuff bcftools ensembl-vep cmake ucsc-gtftogenepred ngs_utils
conda activate pvac

######
# VEP
mkdir vep_data && vep_install -a cf -s homo_sapiens -y GRCh38 -c vep_data/GRCh38
vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream 

pip install pvactools
pvacseq install_vep_plugin /home/563/vs2870/.vep/Plugins

########
# IEDB
tar -zxvf IEDB_MHC_I-2.19.1.tar.gz
# need python2 to run configure :(
conda deactivate 
conda create -n py2 python=2.7
conda activate py2
cd ~/bin
ln -s /g/data3/gx8/extras/vlad/miniconda/envs/py2/bin/python2.7
cd -
cd mhc_i
./configure

sed -i '/import pkg_resources/d' method/netmhc-4.0-executable/netmhc_4_0_executable/__init__.py
sed -i '/import pkg_resources/d' method/netmhcpan-2.8-executable/netmhcpan_2_8_executable/__init__.py
sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhccons-1.1-executable/netmhccons_1_1_executable/bin/pseudofind
sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhc-3.4-executable/netmhc_3_4_executable/netMHC

tar -xzvf IEDB_MHC_II-2.17.5.tar.gz
cd mhc_ii
./configure.py

###########
# INTEGRATE-Neo (for fusionBedpeAnnotator and fusionBedpeSubsetter scripts)
git clone https://github.com/ChrisMaherLab/INTEGRATE-Neo
cd INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1
module load GCC/6.4.0-2.28
./install.sh -o $CONDA_PREFIX/bin

wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.86.gtf.gz Homo_sapiens.GRCh38.86.genePred
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

#######################
## Loading on Raijin ##
cd /g/data3/gx8/projects/Saveliev_pVACtools
source load.sh

###########
## Usage ##
snakemake -p -s Snakefile pvacseq directory=pvac \
    --config \
    sample=diploid_tumor \
    dna_bcbio=/g/data3/gx8/projects/Saveliev_pVACtools/diploid/bcbio_hg38/final \
    rna_bcbio=/g/data/gx8/data/pVAC/GRCh37_wts_samples/final \
    rna_sample=Unknown_B_RNA
"""


import os
from os.path import join, abspath, dirname, isfile, basename, splitext
from ngs_utils.bcbio import BcbioProject
from ngs_utils.utils import flatten
from ngs_utils.logger import critical, info, debug, warn, err
from ngs_utils.file_utils import verify_file

shell.executable(os.environ.get('SHELL', 'bash'))
shell.prefix("")

IEDB_DIR = '/g/data3/gx8/projects/Saveliev_pVACtools'

out_dir = join('pvacseq', SNAME)
work_dir = join('pvacseq', SNAME, 'work')

DNA_BCBIO   = config['dna_bcbio']
DNA_SNAME   = config.get('dna_sname', SNAME)
RNA_PROJECT = config.get('rna_project')
RNA_SNAME   = config.get('rna_sample')


################################################
### hg38 DNA Project. Required for HLA typing.

dna_run = BcbioProject(DNA_BCBIO, include_samples=[DNA_SNAME])
assert dna_run.genome_build == 'hg38'
if len(dna_run.batch_by_name) == 0:
    critical(f'Error: could not find a sample with the name {DNA_SNAME}. Check yaml file for available options: {dna_run.bcbio_yaml_fpath}')
# Batch objects index by tumor sample names
batches = [b for b in dna_run.batch_by_name.values() if not b.is_germline() and b.tumor and b.normal]
assert len(batches) == 1
batch = batches[0]


rule all:
    input:
        pvacseq  = dynamic(join(out_dir, 'pvacseq_results' , 'MHC_Class_{mhc_class}', SNAME + '.final.tsv')),
        pvacfuse = dynamic(join(out_dir, 'pvacfuse_results', 'MHC_Class_{mhc_class}', SNAME + '.final.tsv'))
    output:
        pvacseq  = join(out_dir, SNAME + '_pVACseq_MHC{mhc_class}.tsv'),
        pvacfuse = join(out_dir, SNAME + '_pVACfuse_MHC{mhc_class}.tsv')
    shell:
        'ln -s {input.pvacseq} {output.pvacseq} && ln -s {input.pvacfuse} {output.pvacfuse}'
    

EPITOPE_LENGTHS = '8,9,10,11'
PREDICTORS = 'NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign'
IEDB_DIR = '/g/data3/gx8/projects/Saveliev_pVACtools'


def _pvac_cmdl(tool, input, sample, hla_types, output_dir, other_params=''):
    return (f'{tool} run {input} {sample} "$(cat {hla_types})" {PREDICTORS} {output_dir} ' \
            f'-e {EPITOPE_LENGTHS} {PREDICTORS} -t --top-score-metric=lowest --iedb-install-directory {IEDB_DIR} ' \
            f'--net-chop-method cterm --netmhc-stab --exclude-NAs {other_params}')

rule pvacseq:
    input:  vcf        = join(work_dir, 'somatic.PON.VEP.SELECT.vcf'),
            hla_types  = join(work_dir, 'hla_line.txt'),
    output: out_dir    = directory(join(out_dir, 'pvacseq_results')),
            out_files  = dynamic(join(out_dir, 'pvacseq_results', 'MHC_Class_{mhc_class}', SNAME + '.final.tsv'))
    params: trna_vaf = 10,
    shell:  _pvac_cmdl('pvacseq', '{input.vcf}', SNAME, '{input.hla_types}', output.out_dir, '{params.trna_vaf}')

rule pvacfuse:
    input:  bedpe     = join(work_dir, 'pizzly', rna_sname + '.anno.subset.bedpe'),
            hla_types = join(work_dir, 'hla_line.txt'),
    output: out_dir   = directory(join(out_dir, 'pvacfuse_results')),
            out_files  = dynamic(join(out_dir, 'pvacfuse_results', 'MHC_Class_{mhc_class}', SNAME + '.final.tsv'))
    params: sample = sname,
    shell:  _pvac_cmdl('pvacfuse', '{input.bedpe}', SNAME, '{input.hla_types}', output.out_dir)


##########################
### for both pVAC tools

t_optitype = join(batch.tumor.dirpath, batch.tumor.name + '-hla-optitype.csv')

rule prep_hla:
    """ Output file contains a single line of comma-separated HLA alleles
    """
    input:  t_optitype
    output: join(work_dir, 'hla_line.txt')
    shell:  "cat {input} | grep -v ^sample | tr ',' '\t' | cut -f3 | tr ';' '\n' | sort -u | tr '\n' ',' | head -c -1 > {output}"


################
### for pVACseq
t_vcf = verify_file(join(dna_run.date_dir, f'{batch.name}-ensemble-annotated.vcf.gz'))

rule vcf_pon:
    input:  t_vcf
    output: join(work_dir, 'somatic.PON.vcf')
    params: genome = run.genome_build
    shell:  'pon_anno -g {params.genome} {input} -h1 | bcftools view -f.,PASS -o {output}'
    
rule vcf_vep:
    input:  join(work_dir, 'somatic.PON.vcf')
    output: join(work_dir, 'somatic.PON.VEP.vcf')
    params: assembly = 'GRCh38' if run.genome_build == 'hg38' else 'GRCh37'
    shell:  '''
    vep --input_file {input} --format vcf --output_file {output} --vcf \
    --pick --symbol --terms SO --plugin Downstream --plugin Wildtype \
    --cache --dir_cache ../vep_data/{params.assembly} --assembly {params.assembly} --offline
    '''

rule vcf_select:
    input:  join(work_dir, 'somatic.PON.VEP.vcf')
    output: join(work_dir, 'somatic.PON.VEP.SELECT.vcf')
    params: sample = sname
    shell:  'bcftools view -s {params.sample} {input} > {output}'

t_bam = verify_file(batch.tumor.bam, silent=True)
if t_bam:
    pass
    # Check coverage?

if RNA_SNAME:
    rna_run = BcbioProject(rna_bcbio, include_samples=[RNA_SNAME])
    rna_sample = rna_run.samples[0]
    rna_bam = rna_sample.bam
    rna_tpkm = None
    rna_pizzly_prefix = join(rna_sample.dirpath, 'pizzly', RNA_SNAME)

    rule rna_coverage:
        input:  rna_bam
        output: join(work_dir, '')
        shell:  ''

    rule rna_expression:
        input:  rna_bam
        output: join(work_dir, '')
        shell:  ''

    #################
    ### for pVACfuse

    rule rna_pizzly_to_bedpe:
        input:  rna_pizzly_prefix + '-flat-filtered.tsv', 
                rna_pizzly_prefix + '.json', 
                rna_pizzly_prefix + '.fusions.fasta',
        params: prefix = rna_pizzly_prefix
        output: join(work_dir, 'pizzly', SNAME + '.bedpe')
        shell:  'pizzly_to_bedpe.py {params.prefix} -o {output}'

    rule rna_anno_bedpe:
        input:  join(work_dir, 'pizzly', SNAME + '.bedpe')
        output: join(work_dir, 'pizzly', SNAME + '.anno.bedpe')
        shell:  """
    fusionBedpeAnnotator \
    --reference-file Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --gene-annotation-file Homo_sapiens.GRCh37.75.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file CCR180081_MH18T002P053_RNA-flat-filtered.bedpe \
    --output-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe
    """

    rule rna_subset_bedpe:
        input:  join(work_dir, 'pizzly', SNAME + '.anno.bedpe')
        output: join(work_dir, 'pizzly', SNAME + '.anno.subset.bedpe')
        shell: """
    fusionBedpeSubsetter \
    --input-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe \
    --rule-file INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt \
    --output-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe
    """





        