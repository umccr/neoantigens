# INTEGRATE-Neo

INTEGRATE-Neo is a gene fusion neoantigen discovering tool using next-generation sequencing data. It is written in C++ and Python.

  - Python
  - Perl
  - awk
  - GCC

If not, please install these languages or tools. You may also need to install some prerequisite tools:

  - [BWA](https://sourceforge.net/projects/bio-bwa)
  - [HLAminer v1.3](http://www.bcgsc.ca/platform/bioinfo/software/hlaminer)
  - [NetMHC v4.0](http://www.cbs.dtu.dk/services/NetMHC/output.php)

HLAminer and NetMHC are also included in the vendor directory here. 

To compile the C++ part of this pipeline, you may need to install [CMAKE](https://cmake.org/)

### Installation

Download INTEGRATE-Neo at https://github.com/ChrisMaherLab/INTEGRATE-Neo.

Run the installation script:

```sh
$ cd INTEGRATE-Neo-V-1.2.0
$ chmod +x install.sh
$ ./install.sh -o /opt/bin/
```

Note that you can choose wherever you like to install the software. It can be different from "/opt/bin/". 

Now you have installed:

  - integrate-neo

together with the modules of integrate-neo that can be used as standalone tools:
  - fusionBedpeAnnotator
  - fusionBedpeSubsetter
  - runHLAminer
  - HLAminerToTsv
  - runAddNetMHC4Result
  - runNetMHC4WithSMCRNABedpe

A setup.ini and a rule.txt file are also at your destination directory now. If you don't like them to be there, copy them to the place you like. But remember to use the --setup-file and --rule-file options to run integrate-neo if you moved them.

### setup

Remember to edit the setup.ini file before your first running the pipeline. The one in the installation packages are using example paths like "/SOME/PATH/...".

For the HLAminer reference HLA alleles, i.e. HLA_ABC_CDS.fasta, remember to index it with bwa before the first run.

### input

If you type the following (or python ./integrate-neo.py --help): 

```sh
$ ./integrate-neo.py
```
you can see the 14 parameters and explanations. 

The following are the required options:

        -1/--fastq1       
        -2/--fastq2       
        -f/--fusion-bedpe 
        -r/--reference    
        -g/--gene-model   

The --fastq[1/2] and --reference options are clear enough, the FASTQ and FASTA formats for sequencing reads and human reference genome. 

The --fusion-bedpe option requires a BEDPE format for gene fusions. This BEDPE format follows the standardized format provided by The ICGC-TCGA DREAM Somatic Mutation Calling - RNA Challenge ([SMC-RNA](http://dreamchallenges.org/)).

The --gene-model option requires a gene annotation genePhred file.
  
Download the gtf file from Ensembl:

GRCh37, e.g., v75: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

GRCh38, e.g., v86: ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz

and run the following command for v75:

```sh
$ gunzip Homo_sapiens.GRCh37.75.gtf.gz
$ ./gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.genePred
```

for v86:

```sh
$ gunzip Homo_sapiens.GRCh38.86.gtf.gz.
$ ./gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.86.gtf Homo_sapiens.GRCh38.86.genePred
```

FASTA files can also be downloaded at Ensembl:

v75: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz

v86: ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

### output

The output is in BEDPE format, the first 11 columns follows the SMC-RNA format. columns 12-19 are:
 
 - Epitope sequence
 - Epitope Affinity (nanoMolar)	
 - HLA allele	
 - HLA category	
 - HLA score	
 - HLA e-value	
 - HLA confidence

### Important

The chromosome names in the reference genome, the gene models, and the fusions should be consistent. 

### Examples

Examples are provided for you to test the code.

### Enjoy!

### Release notes:

12-23-2016: INTEGRATE-Neo v 1.2.0

updated BedpeAnnotator to v 0.2.0, which includes a new column for transcript Ids, a new column for lengths of nucleotides in the coding regions at 5p transcripts, a new column for whether the peptides are in-frame, and a new column for whether the fusion transcript follows canonical dinucleotides. 

01-17-2017: INTEGRATE-Neo v 1.2.1

updated BedpeAnnotator to v 0.2.1, which includes a bug fixing.
  


