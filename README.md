# Neoantigen identification

A cancer cell can be targeted by the immune system based on novel mutations ("Variant Antigens"). DNA harbouring such mutations translate into neoantigen peptides. From NGS somatic variant calling data, we can attempt to reconstruct epitopes  sequences of such neoantigen peptides, or more specifically, epitopes, and thus produce personalised DNA-based cancer vaccines. 


# pVACtools

[GitHub](https://github.com/griffithlab/pVACtools)

[Docs](https://pvactools.readthedocs.io/en/latest/pvacseq/prerequisites.html)

[Publication](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0264-5)

pVACtools is a set of 3 pipelines:

* pVACseq - Identifying and prioritizing neoantigens from a list of tumor mutations.
* pVACfuse -	Detecting neoantigens resulting from gene fusions.
* pVACvector - Aid in the construction of DNA-based cancer vaccines.


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

Working on raijin at:

```
/g/data3/gx8/projects/Saveliev_pVACtools
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

Getting error: `FileNotFoundError: [Errno 2] No such file or directory: 'tests/test_data/pvacseq/genes.fpkm_tracking'` Okay, the input file `additional_input_file_list.yaml` has hardcoded paths. Fixing it and rerunning:

```
sed -i s,tests/test_data/pvacseq,pvacseq_example_data,  pvacseq_example_data/additional_input_file_list.yaml
rm -rf pvacseq_example_output  # clean previous output directory
```

Success now. 

## Running offline

However, same command fails on Rajijn worker nodes, because it needs the internet connection.
Maybe install IEDB locally? Because [see at](https://pvactools.readthedocs.io/en/latest/pvacseq/run.html):

> Using a local IEDB installation is strongly recommended for larger datasets or when the making predictions for many alleles, epitope lengths, or prediction algorithms. More information on how to install IEDB locally can be found on the Installation page.

Downloading from [http://tools.iedb.org/mhcii/download](http://tools.iedb.org/mhcii/download/):
clicking `MHC Class I`, then `To download the tools in tar.gz format: Agree and Download`, and then the same for `MHC Class II`. Uncompressing and installing:

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

Now trying to rerun pVACseq with `--iedb-install-directory`:

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

Forgetting the idea to run on worker nodes. As we see further, the login node is sufficient to run the tool on full data. Also keeping in mind an option to fall back to the Dockerized version to run on AWS (see below).


## Running on NeverResponder

### Preparing minimal inputs (a VCF)

The minimal required input is a VCF with VEP-annotated tumor mutations. From umccrised bcbio output, taking post-PoN ensemble somatic VCF, annotating with VEP and bcftool'ing to extract the tumor sample. 

Trying to run on the NeverResponder sample (`/data/cephfs/punim0010/data/Results/Patients/2018-01-17`), as it has both mutation and expression data.

Copying the VCF back from spartan:

```
scp -r spa:/data/cephfs/punim0010/data/Results/Patients/2018-01-17/umccrised/diploid__diploid_tumor/small_variants/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.vcf.gz .
```

Installing VEP:

```
conda install -y -c bioconda ensembl-vep
mkdir vep_data
vep_install -a cf -s homo_sapiens -y GRCh37 -c vep_data/GRCh37
```

Installing [required VEP plugins](https://pvactools.readthedocs.io/en/latest/pvacseq/prerequisites.html):

```
# Installs public Downstream.pm plugin:
vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream  
# Installs pVAC's Wildtype.pm plugin:
pvacseq install_vep_plugin /home/563/vs2870/.vep/Plugins
```

Annotating the VCF:

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

Getting bcftools and extracting the tumor sample:

```
conda install -y -c bioconda bcftools
bcftools view -s diploid_tumor data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.vcf > data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf
```


### Running with VCF

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

Works just nicely on the login node, so hopefully we are not gonna need to set it up for offline work.

### Exploring results

The tool created 2 subfolders in the output directory:

```
/g/data3/gx8/projects/Saveliev_pVACtools/diploid_tumor_output/MHC_Class_I
/g/data3/gx8/projects/Saveliev_pVACtools/diploid_tumor_output/MHC_Class_II
```

Each folder, contains the following files:

* `diploid_tumor.tsv`:                 an intermediate file with variant, transcript, coverage, vaf, and expression information parsed from the input files.
* `diploid_tumor.combined.parsed.tsv`:  above + with binding scores from IEDB added.
* `diploid_tumor.filtered.binding.tsv`: above + filtering by binding threshold.
* `diploid_tumor.stab.tsv`:             above + stability predictions added (optional).
* `diploid_tumor.final.tsv`: 	         The final output file after all filtering and optional steps.


For `MHC_Class_I`, that's one final epitope in `diploid_tumor.final.tsv`:

```
Chromosome  Start     Stop      Reference  Variant  Transcript       Ensembl Gene ID  Variant Type  Mutation  Protein Position  Gene Name  HLA Allele   Peptide Length  Sub-peptide Position  Mutation Position  MT Epitope Seq  WT Epitope Seq  Best MT Score Method  Best MT Score  Corresponding WT Score  Corresponding Fold Change  Tumor DNA Depth  Tumor DNA VAF  Tumor RNA Depth  Tumor RNA VAF  Normal Depth  Normal VAF  Gene Expression  Transcript Expression  Median MT Score  Median WT Score  Median Fold Change  NetMHC WT Score  NetMHC MT Score  PickPocket WT Score  PickPocket MT Score  Best Cleavage Position  Best Cleavage Score  Cleavage Sites                                          Predicted Stability  Half Life  Stability Rank  NetMHCstab allele
14          92470780  92470781  C          T        ENST00000267622  ENSG00000100815  missense      C/Y       1180              TRIP11     HLA-G*01:09  9               10                    2                  KYQTLLAVL       KCQTLLAVL       PickPocket            384.090444123  2216.49442399           5.771                      NA               NA             NA               NA             NA            NA          NA               NA                     384.090444123    2216.49442399    5.771               NA               NA               2216.49442399        384.090444123        9                       0.955652             2:0.804808,5:0.926654,6:0.896452,8:0.599339,9:0.955652  0.058                0.24       0.30            HLA-A24:03
```

`MHC_Class_II` has 106 final epitopes reported, based on 20 mutations in 20 different genes.

If we run with `-t` option, only 1 epitope per variant will be reported, and the output will also contain an intermediate file:

* `diploid_tumor.top.tsv`:  above + picking the top epitope for each variant (optional).

## Adding expression and coverage

Coverage and expression data can be added to the pVACseq processing by providing bam-readcount and/or Cufflinks output files as additional input files. These additional input files must be 
provided as a [yaml file in the following structure](https://pvactools.readthedocs.io/en/latest/pvacseq/prerequisites.html#optional-preprocessing):

```
gene_expn_file:              <genes.fpkm_tracking file from Cufflinks>
transcript_expn_file:        <isoforms.fpkm_tracking file from Cufflinks>
normal_snvs_coverage_file:   <bam-readcount output file for normal BAM and snvs>
normal_indels_coverage_file: <bam-readcount output file for normal BAM and indels>
tdna_snvs_coverage_file:     <bam-readcount output file for tumor DNA BAM and snvs>
tdna_indels_coverage_file:   <bam-readcount output file for tumor DNA BAM and indels>
trna_snvs_coverage_file:     <bam-readcount output file for tumor RNA BAM and snvs>
trna_indels_coverage_file:   <bam-readcount output file for tumor RNA BAM and indels>
```

### Preparing coverage

Splitting the variants into SNPs and indels, and converting to site-list files for bam-readcount:

```
bcftools view -v snps diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf | bcftools query -f "%CHROM\t%POS\n" | awk '{OFS="\t"}{ print $1,$2,$2 }' > diploid/tumor_snps_sites.txt
bcftools view -v indels diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf | bcftools query -f "%CHROM \t %POS \t %LEN \n" | awk '{OFS="\t"}{ print $1,$2,$2+$3-1 }' > diploid/tumor_indels_sites.txt
```

Installing bam-readcount with conda into the `pvac` environment causes a conflict with `perl-storable`, thus creating a new environment with `conda create -n bam-readcount bam-readcount`.

BAM files for coverage:

Tumor DNA:               `/g/data3/gx8/projects/Saveliev_Diploid/diploid_tumor-ready.bam`
Normal DNA (spartan):    `/data/cephfs/punim0010/data/Results/Patients/2018-01-17/final/diploid_blood/diploid_blood-ready.bam`
RNA to genome (spartan): `/data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam`

Running bam-readcount for the NeverResponder (`-w 1` to suppress repeated warnings, `-i` for indels):

```
# tumor sample
bam-readcount -f /g/data3/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa           -w 1 -l diploid/tumor_snps_sites.txt   diploid_tumor-ready.bam    > diploid/tumor_snps.readcount
bam-readcount -f /g/data3/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa           -w 1 -l diploid/tumor_indels_sites.txt diploid_tumor-ready.bam -i > diploid/tumor_indels.readcount

# on spartan 
cd /data/cephfs/punim0010/extras/vlad/tmp
# blood sample
bam-readcount -f /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -w 1 -l tumor_snps_sites.txt           diploid_blood-ready.bam -i > normal_snps.readcount
bam-readcount -f /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -w 1 -l tumor_indels_sites.txt         diploid_blood-ready.bam -i > normal_indels.readcount
# rna sample
bam-readcount -f /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -w 1 -l tumor_snps_sites.txt   /data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam    > rna_snps.readcount
bam-readcount -f /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -w 1 -l tumor_indels_sites.txt /data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam -i > rna_indels.readcount
```

### Running with coverage

Creating yaml file `data/diploid/additional_input_file_list.coverage.yaml`:

```
normal_snvs_coverage_file:   /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/normal_snps.readcount
normal_indels_coverage_file: /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/normal_indels.readcount
tdna_snvs_coverage_file:     /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/tumor_snps.readcount
tdna_indels_coverage_file:   /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/tumor_indels.readcount
trna_snvs_coverage_file:     /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/rna_snps.readcount
trna_indels_coverage_file:   /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/rna_indels.readcount 
```

Executing pvacseq into `diploid_tumor_output_coverage `:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
diploid_tumor_output_coverage \
-e 9,10 \
--tdna-vaf 20 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.coverage.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools
```

A new file is added to output folder:

* `diploid_tumor.coverage.tsv`:  above + filtering on coverage, VAF, and expression values.

And pVACseq filtered out all variants now in `MHC_Class_II`, and reduced `MHC_Class_I` epitopes count from 106 down to 46 (20->10 snps). That's because of set AF threshold for tumor in DNA (`tdna-vaf 20`) and the default RNA threashold of 40%. E.g. the variant with DNA AF=15% was removed:

```
chromosome_name  start     stop      reference  variant  gene_name  amino_acid_change  variant_type  protein_position  transcript_expression  gene_expression  normal_depth  normal_vaf  tdna_depth  tdna_vaf            trna_depth  trna_vaf
14               92470780  92470781  C          T        TRIP11     C/Y                missense      1180              NA                     NA               70            0.0         113         15.929202130159103  138         34.782606
```

So rerunning to set the threshold to 10%:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
diploid_tumor_output_coverage_af10 \
-e 9,10 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.coverage.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--tdna-vaf 10 \
--trna-vaf 10
```

Now getting 1 epitope in `MHC_Class_I` back, and 73 epitopes (14 variants) in `MHC_Class_II`. It filtered out variants with normal AF>2%, for example a variant with tumor AF 18%, and normal AF 8% was filtered out. I think since we're doing all sensible somatic filtering before passing the VCF into pVACtools, we don't need to do extra filtering here. So we will just remove the frequencies and depths from the inputs (except for those for RNA).  

However, one note: if the read coverage is 0, bam-readcount won't report such SNP at all, and pVACseq will use it as NA, and will keep such SNP as passed. For example, consider the following 2 SNPs. The second one has 2 read coverage, and it's filtered out correctly. The first one has 0 read coverage, and reported as NA, so the epitopes that this SNP supports will pass the filters:

```
chromosome_name  start      stop       reference  variant  gene_name  trna_depth  trna_vaf
1                47691408   47691409   A          C        TAL1       NA          NA
2                179558398  179558399  C          G        TTN        2           0.0
```

Wee don't want such SNPs with NA to be reported, as zero RNAseq coverage means no expression, and no epitopes to target by the immune system. So we should post-filter those downstream, e.g. with:

```
awk -F"\t" '{OFS="\t"} { if ($24=="Tumor RNA Depth" || $24!="NA") print; } ' diploid_tumor.final.tsv 
```


### Expression data

RNA coverage above in principle should be sufficient for expression-based filtering, however pVACseq can make use of fpkm count files from cufflinks, so we will test it with them as well:

Prerequisites (on spartan):

```
conda install -y -c bioconda cufflinks
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools
```

Running against the DNA-aligned BAM:

```
cufflinks -p 30 -o cufflinks /data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam
```

It wouldn't work - need transcript IDs. Passing the GTF?

```
cufflinks -p 30 -o cufflinks_gtf /data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam --GTF /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf
```

Looks good. Creating yaml file `data/diploid/data/diploid/additional_input_file_list.expression.gtf.yaml `:

```
gene_expn_file:            /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/genes.fpkm_tracking
transcript_expn_file:      /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/isoforms.fpkm_tracking
```

Running:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
diploid_tumor_output_expression \
-e 9,10 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10
```

`MHC_Class_II` has 84 epitopes in 14 variants (down from 106 and 20). `MHC_Class_I` has the same one variant passed with seemingly correctly annotated gene and transcript expression:

```
cut -f28-29 diploid_tumor.final.tsv
Gene Expression  Transcript Expression
16.2014          4.10986
```


### Expression + RNA coverage

Inputs file `data/diploid/data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml`:

```
gene_expn_file:            /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/genes.fpkm_tracking
transcript_expn_file:      /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/isoforms.fpkm_tracking
trna_snvs_coverage_file:   /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/rna_snps.readcount
trna_indels_coverage_file: /g/data3/gx8/projects/Saveliev_pVACtools/data/diploid/rna_indels.readcount
```

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
diploid_tumor_output_rnacov_expression \
-e 9,10 \
--net-chop-method cterm \
--netmhc-stab \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10
```

`MHC_Class_II` has now 64 epitopes in 12 variants (down from 106 and 20 originally), the rest are filtered out based either on RNA depth/AF or expression. `MHC_Class_I` has the same one variant passed with seemingly correctly annotated gene and transcript expression:

```
cut -f24-25,28-29 diploid_tumor.final.tsv | tsv
Tumor RNA Depth  Tumor RNA VAF      Gene Expression  Transcript Expression
138.0            34.78260617517346  16.2014          4.10986
```

Running with `-t` to get the top results might even reduce the output to a human readable file.


## Docker and CWL
https://pvactools.readthedocs.io/en/latest/install.html#iedb-install
A Docker container for pVACtools is available on DockerHub using the mgibio/pvactools repo:
https://hub.docker.com/r/mgibio/pvactools

CWL tool wrappers for pVACseq, pVACfuse, and pVACvector can be downloaded using the pvactools download_cwls command. The CWLs do not support the –iedb-install-directory or –additional-input-file-list options. Download CWL tool wrappers:

```
pvactools download_cwls cwls
```



















