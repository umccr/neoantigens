# 🏺 Neoantigen identification with pVACseq <!-- omit in toc --> 

Anti-tumor T cells recognize tumor somatic mutations, translated as single amino acid substitutions, as "neoantigens" (upon processing and presentation in the context of MHC molecules). From NGS somatic variant calling data, we can attempt to reconstruct peptide sequences of such neoantigens, and thus produce personalised cancer vaccines.

- [Installation](#installation)
  - [Docker and CWL](#docker-and-cwl)
- [Running on example data](#running-on-example-data)
- [Running offline](#running-offline)
- [Preparing VCF](#preparing-vcf)
- [Preparing HLA types](#preparing-hla-types)
- [Minimal run](#minimal-run)
- [Adding coverage](#adding-coverage)
  - [Preparing coverage](#preparing-coverage)
  - [Running with coverage](#running-with-coverage)
- [Expression data](#expression-data)
  - [Expression from other tools](#expression-from-other-tools)
  - [Expression + RNA coverage](#expression--rna-coverage)
- [Other options](#other-options)
  - [Prediction algorithms](#prediction-algorithms)
  - [Scoring metric](#scoring-metric)
  - [NetChop prediction](#netchop-prediction)
  - [NetMHCstabpan prediction](#netmhcstabpan-prediction)
  - [Downstream limit for frameshifts](#downstream-limit-for-frameshifts)
  - [Top epitope](#top-epitope)
- [HLA typing](#hla-typing)
- [Command line for production](#command-line-for-production)
- [Downstream analysis](#downstream-analysis)
  - [pVACvector](#pvacvector)
  - [pVACviz](#pvacviz)
- [Results validity](#results-validity)

pVACseq is a workflow that identifies and shortlist candidate neoantigen peptides from tumor mutational repertoire. It can use both DNA and RNA sequencing data. Predicted peptides can be used in a personalized vaccine after immunological screening. The tool offers the functionality to compare and differentiate the epitopes found in normal cells against the neoepitopes specifically present in tumor cells for use in personalized cancer vaccines.

[GitHub](https://github.com/griffithlab/pVACtools)

[Docs](https://pvactools.readthedocs.io/en/latest/pvacseq/prerequisites.html)

[Publication](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0264-5)

pVACtools is a set of 3 pipelines:

* **pVACseq**: Identifies and prioritizes neoantigens from a list of tumor mutations. Takes somatic VCF and HLA types, and 1. performs epitope prediction, 2. integrates sequencing-based information, 3. filters neoantigen candidates.
* **pVACfuse**:	Detects neoantigens resulting from gene fusions.
* **pVACvector**: Aids in the construction of DNA-based cancer vaccines.

![Structure](pVACtools_main-figure_v2e.png)


## Installation

Working on Raijin.

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


### Docker and CWL

A [Docker container](https://pvactools.readthedocs.io/en/latest/install.html#docker-and-cwl) for pVACtools is available on DockerHub using the [mgibio/pvactools repo](https://hub.docker.com/r/mgibio/pvactools).

CWL tool wrappers for pVACseq, pVACfuse, and pVACvector can be downloaded using the `download_cwls` command. The CWLs do not support the `--iedb-install-directory` or `--additional-input-file-list` options. Download CWL tool wrappers:

```
pvactools download_cwls cwls_dir
```


## Running on example data

Following the [docs](https://pvactools.readthedocs.io/en/latest/pvacseq/getting_started.html), download the example data (will be saved into `./pvacseq_example_data`):

```
pvacseq download_example_data .
```

And running:

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

Processes IEDB correctly, however still fails because needs the internet for NetChop and NetMHCstabpan predictions:

```
requests.exceptions.ConnectionError: HTTPConnectionPool(host='www.cbs.dtu.dk', port=80)
```

We can turn off running the NetChop and NetMHCstabpan predictions by removing the options `--net-chop-method` and `--netmhc-stab`:

```
pvacseq run \
pvacseq_example_data/input.vcf \
Test \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NetMHC PickPocket NNalign \
pvacseq_example_output_localdb_nonetchop_nostab \
-e 9,10 \
-i pvacseq_example_data/additional_input_file_list.yaml \
--tdna-vaf 20 \
--top-score-metric=lowest \
-d full \
--keep-tmp-files \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools
```

Works offline now!

As we see further, `NetChop` and `NetMHCstabpan ` just add some annotations and don't seem to do any filtering (however those can be done manually), so are quite optional. On the other side, as we see further, the login node is sufficient to run the tool on full data. Also keeping in mind an option to fall back to the Dockerized version to run on AWS.


## Preparing VCF

The tool's minimal required inputs are:

- VCF of somatic mutations, annotated with amino acid changes and transcript sequences by VEP.
- HLA haplotypes of the patient.

The VCF file should contain VEP-annotated tumor mutations. From umccrised bcbio output, we are taking post-PoN ensemble somatic VCF, then annotating with VEP and bcftool'ing to extract the tumor sample.

We are using NeverResponder sample (`/data/cephfs/punim0010/data/Results/Patients/2018-01-17`), as it has both mutation and expression data.

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

## Preparing HLA types

We need to know HLA alleles to properly run pVACseq. 

The human leukocyte antigen (HLA) cluster located on chromosome 6 is one of the most polymorphic regions of the human genome and encodes for several genes involved in functions of the immune system, including HLA classes I and II. Both HLA classes comprise three major loci (HLA-I: A, B, C; HLA-II: DP, DQ, DR), which are co-dominantly expressed. HLA-I/II molecules present intracellular and extracellular peptides, respectively, and interact with other immune cells to induce an adaptive immune response. 

HLA typing can be done with different degrees of resolution, with two-digit and four-digit types distinguishing HLA allele families and distinct HLA protein sequences, respectively. Typing can be done either using specific probing techniques, or purely in-silico from NGS sequencing data through computational means. However, because of the high variability of the HLA loci, the typical read mapping and variant calling-based analysis of NGS data is not suitable to determine the HLA genotype.

[pVACseq methods](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0264-5#CR11) section desribes their appoach to to HLA typing:
- They used clinically genotyped calls whenever they are available,
- For in-silico, they typed on the normal (peripheral blood mononuclear cells), rather than the tumor sample,
- They used 2 tools (HLAminer or by Athlates) and note that they were >85% concordant, but it is helpful to use both algorithms in order to break ties reported by HLAminer.
- Some epitope prediction algorithms, including NetMHC [13](https://scholar.google.com/scholar?hl=en&q=Lundegaard%20C%2C%20Lamberth%20K%2C%20Harndahl%20M%2C%20Buus%20S%2C%20Lund%20O%2C%20Nielsen%20M.%20NetMHC-3.0%3A%20accurate%20web%20accessible%20predictions%20of%20human%2C%20mouse%20and%20monkey%20MHC%20class%20I%20affinities%20for%20peptides%20of%20length%208-11.%20Nucleic%20Acids%20Res.%202008%3B36%28Web%20Server%20issue%29%3AW509%E2%80%93512.), [14](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=12717023), only work with an algorithm-specific subset of HLA alleles, so we are constrained to the set of NetMHC-compatible alleles (e.g. NetMHC v3.4 supports 78 human alleles). On the other hand, such specific epitope prediction software perform slightly better when compared to pan-specific methods such as NetMHCpan in case of well-characterized alleles due to availability of large amounts of training data.

[Bcbio supports](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html?highlight=hla#hla-typing) two HLA callers: `optitype` and `bwakit`. In any case, it first makes use of [bwakit](https://github.com/lh3/bwa/tree/master/bwakit)'s  `bwa-postalt.js` script in order to [extract HLA reads](https://github.com/bcbio/bcbio-nextgen/blob/a3473775db06540c10b5f20ddc2043b8cc99d1f8/bcbio/ngsalign/bwa.py#L70). The different is that the result is passed either to [OptiType's](https://github.com/FRED-2/OptiType) `OptiTypePipeline.py`, or to bwakit's run-HLA. 

According to OptiType (2014) benchmarks, HLAminer (2012) consistently underperformes compared to OptiType, and Athlates (2013) show comparable results.

Bcbio does HLA typing only from hg38 alignments. That's probably because hg38 contains HLA alleles as alternative assembly contigs, which can significantly aid the typing. We will re-analyse the sample against the hg38 build with bcbio with enabled HLA typing stage. Since bcbio can't run both callers at the same time, we will try `optitype`, but keep in mind that we might also run with `bwakit` for control. 

```
# from spartan:
# copy fastq /data/cephfs/punim0010/data/FASTQ/171220_A00130_0036_BH32JNDSXX/PRJ170218A_SFRC01059_T_* /data/cephfs/punim0010/data/FASTQ/171220_A00130_0036_BH32JNDSXX/PRJ170198_SFRC01059_B_*`
# copy samples csv from previous grch37 run on Spartan /data/cephfs/punim0010/data/Results/Patients/2018-01-17/config/diploid.csv
# copy standard cancer workflow template /g/data3/gx8/projects/std_workflow/std_workflow_cancer.yaml and add `hlacaller: optitype` into it
# copy pbs submitter /g/data3/gx8/projects/std_workflow/run.sh
# running bcbio on spartan:
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools/diploid/bcbio_hg38/diploid/work
source ~/load_bcbio.sh
bcbio_nextgen.py -w template std_workflow_cancer.yaml diploid.csv *.fastq.gz
cp run.sh diploid/work
cd diploid/work
qsub run.sh
```

Optitype produced the following alleles (from `diploid_blood-optitype.csv` and `duploid_tumor-optitype.csv`), identical for tumor and blood:

```
locus  alleles                
A      HLA-A*02:01;HLA-A*26:01
B      HLA-B*35:02;HLA-B*18:01
C      HLA-C*04:01;HLA-C*05:01
```

The resulting alleles can be checked against the list of valid alleles for pVACseq:

```
pvacseq valid_alleles > valid_alleles
for A in "HLA-A*02:01" "HLA-A*26:01" "HLA-B*35:02" "HLA-B*18:01" "HLA-C*04:01" "HLA-C*05:01" ; do grep -q -F $A valid_alleles; done
# should exit code zero
```

All 6 alelles are valid.

Most of the tests below we did before we had HLA typing results from bcbio run against the hg38 build, so we used common `HLA-G*01:09`, `HLA-E*01:01`, and `H2-IAb` as input types instead. [Later in this article we tested this sample with proper HLA types as well](#hla-typing).


## Minimal run 

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
* `diploid_tumor.final.tsv`: 	         The final output file after all filtering and optional steps (would be the same as `diploid_tumor.stab.tsv`).

For `MHC_Class_I`, that's one final epitope in `diploid_tumor.final.tsv`:

```
Chromosome  Start     Stop      Reference  Variant  Transcript       Ensembl Gene ID  Variant Type  Mutation  Protein Position  Gene Name  HLA Allele   Peptide Length  Sub-peptide Position  Mutation Position  MT Epitope Seq  WT Epitope Seq  Best MT Score Method  Best MT Score  Corresponding WT Score  Corresponding Fold Change  Tumor DNA Depth  Tumor DNA VAF  Tumor RNA Depth  Tumor RNA VAF  Normal Depth  Normal VAF  Gene Expression  Transcript Expression  Median MT Score  Median WT Score  Median Fold Change  NetMHC WT Score  NetMHC MT Score  PickPocket WT Score  PickPocket MT Score  Best Cleavage Position  Best Cleavage Score  Cleavage Sites                                          Predicted Stability  Half Life  Stability Rank  NetMHCstab allele
14          92470780  92470781  C          T        ENST00000267622  ENSG00000100815  missense      C/Y       1180              TRIP11     HLA-G*01:09  9               10                    2                  KYQTLLAVL       KCQTLLAVL       PickPocket            384.090444123  2216.49442399           5.771                      NA               NA             NA               NA             NA            NA          NA               NA                     384.090444123    2216.49442399    5.771               NA               NA               2216.49442399        384.090444123        9                       0.955652             2:0.804808,5:0.926654,6:0.896452,8:0.599339,9:0.955652  0.058                0.24       0.30            HLA-A24:03
```

`MHC_Class_II` has 106 final epitopes reported, based on 20 mutations in 20 different genes.

The columns are nicely explained [here](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html).


## Adding coverage

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

Now getting 1 epitope in `MHC_Class_I` back, and 73 epitopes (14 variants) in `MHC_Class_II`. It filtered out variants with normal AF>2%, for example a variant with tumor AF 18%, and normal AF 8% was filtered out. I think since we're doing all sensible somatic filtering before passing the VCF into pVACtools, we don't need to do extra filtering here. For example, filtering out normal AF>2% will have trouble with tumor contamination, and want to deal with that before running pVAC. So we will just remove the frequencies and depths from the inputs (except for those for RNA).

However, keeping in mind that we want to use only mutations for building epitopes that consistently appear in a tumor altogether. It's reached ideally on the VAF roughly equal to tumor ploidy. That's why it's good to figure out tumor ploidiy based on typical somatic mutations AFs, and then filter out mutations with tumor VAF below that value.

Another technical note: if the read coverage is 0, bam-readcount won't report such SNP at all, and pVACseq will use it as NA, and will keep such SNP as passed. For example, consider the following 2 SNPs. The second one has 2 read coverage, and it's filtered out correctly. The first one has 0 read coverage, and reported as NA, so the epitopes that this SNP supports will pass the filters:

```
chromosome_name  start      stop       reference  variant  gene_name  trna_depth  trna_vaf
1                47691408   47691409   A          C        TAL1       NA          NA
2                179558398  179558399  C          G        TTN        2           0.0
```

Wee don't want such SNPs with NA to be reported, as zero RNAseq coverage means no expression, and no epitopes to target by the immune system. So we should post-filter those downstream, e.g. with:

```
awk -F"\t" '{OFS="\t"} { if ($24=="Tumor RNA Depth" || $24!="NA") print; } ' diploid_tumor.final.tsv 
```

UPD: can use `--exclude-NAs` option, actually.


## Expression data

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

### Expression from other tools

It's possible to use other bcbio output as expression input too. This is what file formats are actually needed for pVACseq ([from FAQ](https://pvactools.readthedocs.io/en/latest/pvacseq/frequently_asked_questions.html)]:

> For transcript FPKM: a tab-separated file with a tracking_id column containing Ensembl transcript IDs and a FPKM column containing FPKM values.

> For gene FPKM: a tab-separated file with a tracking_id column containing Ensembl gene IDs, a locus column describing the region within the gene, and a FPKM column containing FPKM values. In the pVACseq pipeline the FPKM values will be summed for all loci of a gene. You may also provide already summed FPKM values. In that case you will still need to provide a locus column but the values in that column can be empty.


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

The run took 17 minutes, which is acceptable.


## Other options

We used the following options derived from the example run:

* `HLA-G*01:09,HLA-E*01:01,H2-IAb` (HLA alleles)
* `NetMHC PickPocket NNalign` (Prediction algorithms)
* `-e 9,10` (Length of subpeptides (neoepitopes) to predict. Typical epitope lengths vary between 8-11. Required for Class I prediction algorithms)
* `--net-chop-method cterm` (Use NetChop, and the value is prediction method - `cterm` or `20s`)
* `--netmhc-stab` (Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes)
* `--top-score-metric=lowest` (The ic50 scoring metric to use when filtering epitopes
                        by binding-threshold or minimum fold change. Default: median)
* `-d full` (Cap to limit the downstream sequence length for
                        frameshifts when creating the fasta file. Use 'full'
                        to include the full downstream sequence. Default: 1000)

We will attempt few runs with options derived from those, following the [options docs](https://pvactools.readthedocs.io/en/latest/pvacseq/run.html).


### Prediction algorithms

The following list of the algorithms are available:

```
NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC SMMalign
```

We used `NetMHC PickPocket NNalign`, however we will try to use all of them now:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC SMMalign \
diploid_tumor_output_rnacov_expression_methods \
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

Getting more output now. 

`MHC_Class_I` contains now 5 epitopes, for 5 variants in 5 genes, 3 of which are best predicted with `NetMHCpan`, 1 with `SMM`, and one with `SMMPMBEC`. The one on 14:92470780 that was predicted with `PickPocket` before, is now also best predicted with `NetMHCpan`.

`MHC_Class_II` contains 75 epitopes now, in 13 variants.

Using all methods since it doesn't add much to the computational time (still 26 minutes) and doesn't clutter much the output, and it looks like it gives more reliable prediction overall. However [not all methods are independent](https://pvactools.readthedocs.io/en/latest/pvacseq/frequently_asked_questions.html): `NetMHCcons` is a consensus method between `NetMHC`, `NetMHCpan`, and `PickPocket`, so we can omit those and stick to the following:

```
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign
```


### Scoring metric

Since many algorithms can predict an epitop, the option `--top-score-metric` controls which algorithm's score will be used to report it, and to filter (with `-b` for binding threshold, which is default 500). Two options to select score from all prediction methods:

* `lowest`: the best one (also recorded in the output as `Best MT Score`),
* `median`: median one (`Median MT Score`)

We used `lowest`; now trying with the default `median` (and also setting to use all methods for more representative results):

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC SMMalign \
diploid_tumor_output_rnacov_expression_methods_median \
-e 9,10 \
--top-score-metric=median \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10
```

The output is briefer, with `MHC_Class_II` having 34 epitopes (in 5 variants), and `MHC_Class_I` having 1 epitope. So basically filters applied not to the best but median score drops out some variants. Might wanna stick to the `lowest` as a more safe one, because not sure about the methods weight (especially with collapsing all MHC methods into `NetMHCcons`).


### NetChop prediction 

[NetChop](http://www.cbs.dtu.dk/services/NetChop/) is a server that produces neural network predictions for human proteasome cleavage sites.

If the option `--net-chop-method` is omitted, the NetChop predictions won't run, which saves a couple of minutes of runtime. Turning on NetChop doesn't seem to filter more output, but it produced extra columns as:

```
Best Cleavage Position  Best Cleavage Score  Cleavage Sites
9                       0.955652             2:0.804808,5:0.926654,6:0.896452,8:0.599339,9:0.955652
```

We used `cterm` method. Available are `cterm` or `20s`:

* `cterm` is C-term 3.0 network which is trained on 1260 MHC class I ligands (using only C-terminal cleavage site of the ligands). C-term 3.0 network performs best in predicting the boundaries of CTL epitopes.
* `20s` is 20S network which is trained with in vitro degradation data published in Toes, et al. and Emmerich et al. 

Also note that turning on NetChop requires Internet connectivity. So if we need to run pVACtools offline, we might want to turn it off. Otherwise, sticking to `cterm`.


### NetMHCstabpan prediction

[NetMHCstabpan](http://www.cbs.dtu.dk/services/NetMHCstabpan/) is a server predicts binding stability of peptides to any known MHC molecule using artificial neural networks (ANNs). The method is trained on more than 25,000 quantitative stability data covering 75 different HLA molecules.

Using the `--netmhc-stab` option, like with NetChop, just adds a couple of columns into the output:

```
Predicted Stability  Half Life  Stability Rank  NetMHCstab allele
0.000                0.08       46.00           HLA-B14:01
```

And should be omitted for offline mode.

We might want to consider adding some filtering based on NetChop and NetMHCstabpan results. For example, many epitopes seem to have 0 stability, and we can require some threshold on that.


### Downstream limit for frameshifts

`-d` option limits the downstream sequence length for frameshifts when creating the fasta file. Use `full` to include the full downstream sequence. Default: 1000.

We ran with `d full`. The result is the same with the default value, so sticking with it.


### Top epitope

Running with `-t` reports one top result and makes the output human readable:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC SMMalign \
diploid_tumor_output_rnacov_expression_methods_median_top \
-e 9,10 \
--top-score-metric=median \
-d full \
--keep-tmp-files \
-i data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10 \
-t
```

Getting expected 5 epitopes for 5 variants in `MHC_Class_II`, and 1 for 1 in `MHC_Class_I`.

Also getting 1 additional file in the output folders: `diploid_tumor.top.tsv`, and `diploid_tumor.filtered.coverage.tsv` corresponds to the pre-`-t` output, which is good if we want to go back to other epitopes for a variant, or re-filter (like above)

We can alternatively post-filter already created results with the [top_score_filter command](https://pvactools.readthedocs.io/en/latest/pvacseq/filter_commands.html#top-score-filter):

```
pvacseq top_score_filter -m median MHC_Class_II/diploid_tumor.final.tsv MHC_Class_II/diploid_tumor.final.TOP.tsv
```

The command will reduce the output file to have only 1 epitope per variant. The command is confusing because it uses the same notation (top_score_filter of lowest/median) as the main tool's `-m` option, which selects one prediction across methods. However, this one selects the best across epitops for single variant, reducing 64 records to 12.


## HLA typing

For previous tests, we used generic HLA alleles (taken from example data) because we didn't have true HLA calls at the moment. Now that we [ran bcbio with Optitype against hg38](#hla-types), we can test pVACseq with correct HLA alleles:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01" \
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
diploid_tumor_output_rnacov_expression_hla \
-e 9,10 \
--top-score-metric=lowest \
-i data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10
```

Since all HLA alelels are from MHC class I, we don't have anything for MHC class II. For MHC class 1, we have 88 epitopes predicted, in 31 different variants, from 30 genes.


## Command line for production

Based on all the experiments, this would be the production command line for this sample:

```
pvacseq run \
data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
diploid_tumor \
"HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01" \
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
diploid_tumor_production \
-e 9,10 \
-t \
--top-score-metric=lowest \
-i data/diploid/additional_input_file_list.rna_coverage.expression.gtf.yaml \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10 \
--net-chop-method cterm \
--netmhc-stab
```

For offline runs, we would omit `--net-chop-method cterm` and `--netmhc-stab`

Would keep `-t` to picking only the top epitope for a mutation, keeping in mind that we can go back to intermediate `diploid_tumor.filtered.coverage.tsv`, and refilter with `pvacseq top_score_filter`.

HLA alleles should be repalced for a different sample

Also, with known tumor ploidy, it would be improtant to still filter the mutations with tumor VAF below the ploidy, to make sure the mutations for building epitopes appear in the tumor altogether. In that case, we would add the coverage files back into the inputs yamls and set the `--tdna-vaf` option.


## Downstream analysis

### pVACvector

[pVACvector](https://pvactools.readthedocs.io/en/latest/pvacvector.html) is aids the construction of DNA-based cancer vaccines.

* Input: a pVACseq output.
* Output: ordering that minimizes the effects of junctional epitopes (that may create novel peptides) between the sequences.

It does this by using the core pVACseq services to predict the binding scores for each junctional peptide. 

It also tests junctions with spacer amino acid sequences that may help to reduce reactivity. These spacer amino acid sequences can be “HH”, “HHC”, “HHH”, “HHHD”, “HHHC”, “AAY”, “HHHH”, “HHAA”, “HHL” or “AAL”. The final vaccine ordering is achieved through a simulated annealing procedure that returns a near-optimal solution, when one exists.

The expected output is a fasta file `<sample>_results.fa` with the peptide sequences and best spacers in the optimal order, and a visualization of the above in a PNG file.

Running on `MHC_Class_I` 2 epitopes:

```
pvacvector run \
diploid_tumor_production/MHC_Class_I/diploid_tumor.final.tsv \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
diploid_tumor_production/MHC_Class_I/vector \
-e 9,10 \
--input_vcf data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--keep-tmp-files
```

For some reason it queries the IEDB database again, does it need to? Need to figure out.

The run finished fine. But for some reason didn't produce the expected image - only `.fa` files are found in the output folder. Maybe need to rerun on larger input?


Running on larger `MHC_Class_II`:

```
pvacvector run \
diploid_tumor_production/MHC_Class_II/diploid_tumor.final.tsv \
diploid_tumor \
"HLA-G*01:09,HLA-E*01:01,H2-IAb" \
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
diploid_tumor_production/MHC_Class_II/vector \
-e 9,10 \
--input_vcf data/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--keep-tmp-files
```

Failing with:

```
Running IEDB on Allele HLA-E*01:01 and Epitope Length 9 with Method SMM - Entries 1-2
Traceback (most recent call last):
...
subprocess.CalledProcessError: Command '['python2.7', '/g/data3/gx8/projects/Saveliev_pVACtools/mhc_i/src/predict_binding.py', 'smm', 'HLA-E*01:01', '9', '/g/data3/gx8/projects/Saveliev_pVACtools/diploid_tumor_production/MHC_Class_II/vector/MHC_Class_I/tmp/diploid_tumor_.fa.split_1-2']' returned non-zero exit status -9
```

Hmm, why does it need to predict bindings again? Anyway, trying to rerun the failed command standalone:

```
python2.7 /g/data3/gx8/projects/Saveliev_pVACtools/mhc_i/src/predict_binding.py smm "HLA-E*01:01" 9 /g/data3/gx8/projects/Saveliev_pVACtools/diploid_tumor_production/vector/MHC_Class_I/tmp/diploid_tumor_.fa.split_1-2
```

Fails with OOM. I guess we are gonna need to use the worker nodes this time. Rerunning.

Now it seem to getting stuck on running on of the prediction methods, and running out of walltime:

```
...
Executing MHC Class II predictions
Generating Variant Peptide FASTA and Key Files
Generating Variant Peptide FASTA and Key Files - Entries 1-2
Completed
Processing entries for Allele H2-IAb - Entries 1-2
Running IEDB on Allele H2-IAb with Method NNalign - Entries 1-2
Completed
Running IEDB on Allele H2-IAb with Method NetMHCIIpan - Entries 1-2
=>> PBS: job killed: walltime 36059 exceeded limit 36000
[1]    22797 terminated  pvacvector run
```

Leaving this step out for now.

### pVACviz

Looks like it's unreleased. No information about it on the website, and on GitHub, it seems to be under active development under [staging branch](https://github.com/griffithlab/pVACtools/tree/staging/utils/pvacviz).


## Results validity

In their tests, out of 100-400 missense SNV per sample, they went down to 110-150 epitope candidates with 40-100 in HLA allele of interest, 11-24 after filtering by coverage/expression and eyeballing.

They had 3 samples per patients from different tissue and final epitopes didn't exactly overlap. But they finally chosen 14-18 epitopes to test in binding assays, then going down to 7 to produce vaccines, 3 of which were confirmed immunogenic.

![table](table.png)

