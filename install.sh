####################
### Installation ###
conda env create -p neoantigens --file environment.yml
conda activate neoantigens
pip install -e .

########
## VEP
mkdir vep_data && vep_install -a cf -s homo_sapiens -y GRCh38 -c {VEP_DATA}/GRCh38
vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream 

########
## pVACtools
pip install pvactools
pvacseq install_vep_plugin /home/563/vs2870/.vep/Plugins

########
## For panel of normals- install vcf_stuff

########
## IEDB
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

########
## INTEGRATE-Neo (for fusionBedpeAnnotator and fusionBedpeSubsetter scripts)
cd INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1  # from git clone https://github.com/ChrisMaherLab/INTEGRATE-Neo
module load GCC/6.4.0-2.28  # on raijin
./install.sh -o $CONDA_PREFIX/bin

wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.86.gtf.gz Homo_sapiens.GRCh38.86.genePred
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

########
## pyensembl (to convert pizzly to bedpe)
pip install pyensembl
# In 2 steps: first on loging node to make it download the files:
pyensembl install --release 86 --species human
# when it starts `Reading GTF from`, go into a worker node and run again.

########
## bam-readcount (for pvacseq coverage of rnaseq)
# bam-readcount cannot install into the same conda env, thus we are creating a separate one:
conda create -n bam-readcount -c bioconda -c conda-forge bam-readcount
conda activate bam-readcount
BAM_READCOUNT=$(which bam-readcount)
conda activate neoantigens

#########################
### Loading on Raijin ###
cd /g/data3/gx8/projects/Saveliev_pVACtools
source load.sh

###################
###### Usage ######
snakemake -p -s Snakefile pvacseq --directory pvac \
--config \
dna_sample=diploid_tumor \
dna_bcbio=/g/data3/gx8/projects/Saveliev_pVACtools/diploid/bcbio_hg38/final \
rna_bcbio=/g/data/gx8/data/pVAC/GRCh37_wts_samples/final \
rna_sample=Unknown_B_RNA


## TODO: use the same GTF for:
## - pyensembl
## - fusionBedpeAnnotator
## - cufflinks
## - VEP?
## - bcbio?