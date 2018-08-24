# Neoantigen identification

A cancer cell can be targeted by the immune system based on novel mutations, called Variant Antigens. DNA harbouring such mutations translate into neoantigen peptides. From NGS somatic variant calling data, we can attempt to reconstruct epitopes  sequences of such neoantigen peptides, or more specifically, epitopes, and thus produce personalised DNA-based cancer vaccines. 

# pVACtools
GitHub: https://github.com/griffithlab/pVACtools
Docs: https://pvactools.readthedocs.io/en/latest/pvacseq/prerequisites.html#vep
Publication: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0264-5

A cancer cell can be targeted by the immune system based on novel mutations. We call such mutations Variant Antigens. From NGS somatic variant calls, we can attempt to reconstruct sequences of such neoantigen peptides, or more specifically, epitopes.

pVACtools is a set of 3 tools:
	* pVACseq
	A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a list of tumor mutations.
	* pVACfuse
	A tool for detecting neoantigens resulting from gene fusions.
	* pVACvector
	A tool designed to aid specifically in the construction of DNA-based cancer vaccines.

## Installation (raijin)

```
# source ~/load_conda.sh  # automatically done in ~/.zshrc
conda create -n pvac python=3.5 bcftools
conda activate pvac
pip install pvactools
```

Loading:

```
conda activate pvac
```

## Running on example data

Download the example data (will be saved into `./pvacseq_example_data`)

```
pvacseq download_example_data .
```

Run:

```
pvacseq run \
pvacseq_example_data/input.vcf \
Test \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
pvacseq_example_output2 \
-e 9,10 \
-i pvacseq_example_data/additional_input_file_list.yaml \
--tdna-vaf 20 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files
```

Getting error:
`FileNotFoundError: [Errno 2] No such file or directory: 'tests/test_data/pvacseq/genes.fpkm_tracking'`
Fixing input and rerunning:

```
sed -i s,tests/test_data/pvacseq,pvacseq_example_data,  pvacseq_example_data/additional_input_file_list.yaml
rm -rf pvacseq_example_output  # clean previous output directory
```

Success with that. 

However, it fails on Rajijn worker nodes, because it needs the internet connection :-(
Maybe install IEDB locally? Because see https://pvactools.readthedocs.io/en/latest/pvacseq/run.html :

> Using a local IEDB installation is strongly recommended for larger datasets or when the making predictions for many alleles, epitope lengths, or prediction algorithms. More information on how to install IEDB locally can be found on the Installation page.

Downloading from http://tools.iedb.org/mhcii/download/
Clicking `MHC Class I`, then `To download the tools in tar.gz format: Agree and Download`
Then the same for `MHC Class II`.

```
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
```

Now trying to run with `--iedb-install-directory`:

```
pvacseq run \
pvacseq_example_data/input.vcf \
Test \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
pvacseq_example_output_localdb \
-e 9,10 \
-i pvacseq_example_data/additional_input_file_list.yaml \
--tdna-vaf 20 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools
```

Processes IEDB correctly, however still fails because needs the internet for something else:

```
requests.exceptions.ConnectionError: HTTPConnectionPool(host='www.cbs.dtu.dk', port=80): Max retries exceeded with url: /cgi-bin/webface2.fcgi (Caused by NewConnectionError('<urllib3.connection.HTTPConnection object at 0x7f3a9b931588>: Failed to establish a new connection: [Errno -2] Name or service not known',))
```

Forget that. Will fall back to Dockerized version on AWS if need to run in production (see below).

## Preparing inputs

### VCF

ensemble somatic VCF from bcbio+umccrise, containing 1 sample only. Need to be annotated with VEP. 
Copying the data:

```
scp -r spa:/data/cephfs/punim0010/data/Results/Patients/2018-01-17/umccrised/diploid__diploid_tumor/small_variants/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.vcf.gz .
```

Installing VEP:

```
conda install -y -c bioconda ensembl-vep
mkdir vep_data
vep_install -a cf -s homo_sapiens -y GRCh37 -c vep_data/GRCh37
```

Install VEP plugins:

```
# Installs public Downstream.pm plugin:
vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream  
# Installs pVAC's Wildtype.pm plugin:
pvacseq install_vep_plugin /home/563/vs2870/.vep/Plugins
```

Annotating:

```
vep \
--input_file data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.vcf.gz \
--format vcf \
--pick \
--output_file data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.vcf \
--vcf --symbol --terms SO \
--plugin Downstream \
--plugin Wildtype \
--cache \
--dir_cache vep_data/GRCh37 \
--assembly GRCh37 \
--offline
```

Extracting 1 sample:

```
conda install -y -c bioconda bcftools
bcftools view -s diploid_tumor data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.vcf > data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf
```

### Expression and coverage
Coverage and expression data can be added to the pVACseq processing by providing bam-readcount and/or Cufflinks output files as additional input files. These additional input files must be 
provided as a yaml file in the following structure:

```
gene_expn_file: <genes.fpkm_tracking file from Cufflinks>
transcript_expn_file: <isoforms.fpkm_tracking file from Cufflinks>
normal_snvs_coverage_file: <bam-readcount output file for normal BAM and snvs>
normal_indels_coverage_file: <bam-readcount output file for normal BAM and indels>
tdna_snvs_coverage_file: <bam-readcount output file for tumor DNA BAM and snvs>
tdna_indels_coverage_file: <bam-readcount output file for tumor DNA BAM and indels>
trna_snvs_coverage_file: <bam-readcount output file for tumor RNA BAM and snvs>
trna_indels_coverage_file: <bam-readcount output file for tumor RNA BAM and indels>
```

## Running

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
diploid_tumor_output \
-e 9,10 \
--tdna-vaf 20 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files
# -i data/additional_input_file_list.yaml
```

## Docker and CWL
https://pvactools.readthedocs.io/en/latest/install.html#iedb-install
A Docker container for pVACtools is available on DockerHub using the mgibio/pvactools repo:
https://hub.docker.com/r/mgibio/pvactools

CWL tool wrappers for pVACseq, pVACfuse, and pVACvector can be downloaded using the pvactools download_cwls command. The CWLs do not support the –iedb-install-directory or –additional-input-file-list options. Download CWL tool wrappers:

```
pvactools download_cwls cwls
```
