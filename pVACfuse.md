# Neoepitopes from WTS fusion calls with pVACfuse

- [Neoepitopes from WTS fusion calls with pVACfuse](#neoepitopes-from-wts-fusion-calls-with-pvacfuse)
    - [Running on example data](#running-on-example-data)
    - [Fusions callers](#fusions-callers)
        - [INTEGRATE](#integrate)
        - [Other fusion calers](#other-fusion-calers)
        - [Pizzly](#pizzly)
        - [Oncofuse](#oncofuse)
        - [JAFFA](#jaffa)
        - [Clinker](#clinker)
        - [SQUID](#squid)
        - [INTEGRATE-Neo](#integrate-neo)
            - [fusionBedpeAnnotator](#fusionbedpeannotator)
            - [fusionBedpeSubsetter](#fusionbedpesubsetter)
            - [runNetMHC4WithSMCRNABedpe.py](#runnetmhc4withsmcrnabedpepy)
            - [runAddNetMHC4Result.py](#runaddnetmhc4resultpy)
    - [Running pVACfuse with INTEGRATE-Neo output](#running-pvacfuse-with-integrate-neo-output)
    - [Converting pizzly into INTEGRATE-type bedpe](#converting-pizzly-into-integrate-type-bedpe)
    - [pVACfuse on annotated pizzly](#pvacfuse-on-annotated-pizzly)
    - [INTEGRATE-Neo on pizzly and adding filtering before pVACfuse:](#integrate-neo-on-pizzly-and-adding-filtering-before-pvacfuse)

Working directories:

spartan: `/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion`
raijin: `/g/data3/gx8/projects/Saveliev_pVACtools/fusions`

[pVACfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html) is a part of pVACtools package. It predicts neoantigens produced specifically by gene fusion events.

It requires
* [INTEGRATE-Neo](https://github.com/ChrisMaherLab/INTEGRATE-Neo) predictions,
* ... made out of [INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home) fusions calls.
*  `.genePred` file (we have it in bcbio reference data)

## Running on example data

```
cd /g/data3/gx8/projects/Saveliev_pVACtools/fusions
pvacfuse download_example_data .
pvacfuse run \
pvacfuse_example_data/fusions.bedpe.annot \
Test \
"HLA-A*29:02" \
NetMHC \
pvacfuse_example_data/results \
-e 9 \
--top-score-metric=lowest \
--keep-tmp-files
```

## Fusions callers

Useful wiki: [https://github.com/rdocking/fusebench/wiki/components_and_similar_projects](https://github.com/rdocking/fusebench/wiki/components_and_similar_projects)

### INTEGRATE

[INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home) is very poorly documented and requires tophat alignment. And fails to install.

Prepare annotation file (`spartan`)

```
cd INTEGRATE
PRED=/data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.genePred
cut -f 1-10,12 $PRED > tmp.txt
echo -e "#GRCh37.ensGene.name\tGRCh37.ensGene.chrom\tGRCh37.ensGene.strand\tGRCh37.ensGene.txStart\tGRCh37.ensGene.txEnd\tGRCh37.ensGene.cdsStart\tGRCh37.ensGene.cdsEnd\tGRCh37.ensGene.exonCount\tGRCh37.ensGene.exonStarts\tGRCh37.ensGene.exonEnds\tGRCh37.ensemblToGeneName.value" > annot.txt
cat tmp.txt >> annot.txt
```

Install INTEGRATE:

```
conda install -c bioconda -c conda-forge samtools libdivsufsort cmake bwa
wget https://sourceforge.net/projects/integrate-fusion/files/INTEGRATE.0.2.6.tar.gz/download | tar -xzf -
cd INTEGRATE_0_2_6
mkdir INTEGRATE-build
cd INTEGRATE-build
cmake ../Integrate/ -DCMAKE_BUILD_TYPE=release
make integrate
# doesn't work
```

Doesn't work. Will try figure out what format it needs and make it from kallisto or similar.


### Other fusion calers

pVACfuse expects the input to be a `.bedpe` file as the following. Downloading pVACfuse test data and printing out the exaple fusion file::

(`raijin`)

```
less pvacfuse_example_data/fusions.bedpe.annot
19      -1      39123136        6       -1      46607405        EIF3K>>CYP39A1  1       +       -       2       1       1       QALDENMDLLEGITGFEDSVRKSSIPKNVFLALHEKLYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLV        22      1       ENST00000248342(421)|ENST00000538434(160)|ENST00000545173(421)|ENST00000592558(421)|ENST00000593149(160);ENST00000275016;       1
17      -1      7490589 17      9760739 -1      MPDU1>>GLP2R    3       +       +       4       1       1       LLQAATNYHNGHTGQLSAITVFLLFGGSLARIFTSIQKTPLHAQLHPHELVCFFHPENPGCTGEGRRLLQLLLQEAX,SPGSHQLPQRAHRPALSHHSLPAVWGLPGPNLHFHSENSTARATTSTX,FSRQPPTTTTGTQASSQPSQSSCCLGAPWPESSLPFRKLHCTRNYIHMNLFASFILRTLAVLVKDVVFYNSYSKRPDNENGWMSYLSE 37,36,36        0,1,2   ENST00000250124(618);;ENST00000262441|ENST00000458005|ENST00000574745,ENST00000423172(532);;ENST00000262441|ENST00000458005|ENST00000574745,ENST00000396501(560);ENST00000262441|ENST00000574745;ENST00000458005        1
```

Bcbio RNAseq pipeline supports the following tools:

* kallisto + [pizzly](https://github.com/pmelsted/pizzly)
* STAR + [oncofuse](http://www.unav.es/genetica/oncofuse.html) (GRCh37 only as [STAR does not work with alts in hg38](https://github.com/alexdobin/STAR/issues/39))
* [EricScript](https://sites.google.com/site/bioericscript)

Based on the [discussion](https://github.com/bcbio/bcbio-nextgen/issues/1649), pizzly is the best way to go.

### Pizzly

Examples are at on spartan:

```
/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/pizzly
```

Pizzly output a json file and a flat file.

CCR180109_VPT-EH09_RNA.json (truncated):

```
{
  "genes" : [
    {
      "geneA" : { "id" : "ENSG00000135722", "name" : "FBXL8"},
      "geneB" : { "id" : "ENSG00000204385", "name" : "SLC44A4"},
      "paircount" : 0,
      "splitcount" : 1,
      "transcripts" : [
        {
          "fasta_record": "ENST00000521920_0:255_ENST00000544672_1106:2634",
          "transcriptA": {"id" : "ENST00000521920", "startPos" : 0, "endPos" : 255, "edit" : 0, "strand" : true},
          "transcriptB": {"id" : "ENST00000544672", "startPos" : 1106, "endPos" : 2634, "edit" : 0, "strand" : true},
          "support" : 1,
          "reads" : [0]
        },
        ...
      ],
      "readpairs" : [
        {
          "type" : "SPLIT",
          "read1" : { "name" : "A00130:69:H5VFCDSXX:4:2557:25572:3912", "seq" : "GCATGATGTATTTATTGGCACTTTCTCACAGCCGGGGGAGCTGATGTTGGATGCCCAGAGCACATACTGGGGTTGCCCCGATGTAGCCAGGTACAGAGCAGTCATGGCCCAGTAGGCAATGCAGATGAGGAGGAGGACAAAGGTGACCAGT", "splitpos" : -1, "direction" : "true", "kmerpos" : { "start" : 0, "stop" : 120}},
          "read2" : { "name" : "A00130:69:H5VFCDSXX:4:2557:25572:3912", "seq" : "CCTGTCCCTGAGAGACCGTGCTGCCGCCGCCAGGGTCTGCAGGGCCTGGGCCGCCGCTGCTACCTGCAGCGCCGTGTGGCACGACACAAAAATCAGGGCTGTGGGACAGATGATGTCTACCATGTTCTACCCACTGGTCACCTTTGTCCTC", "splitpos" : 96, "direction" : "false", "kmerpos" : { "start" : 0, "stop" : 65}}
        }
      ]
    },
...
```

The flat file is also filtered on top requiring paircount>=2 and splitcount>=2: 
CCR180109_VPT-EH09_RNA-flat-filtered.tsv (truncated last columns):

```
geneA.name  geneA.id         geneB.name  geneB.id         paircount  splitcount  transcripts.list
TFF1        ENSG00000160182  RPL7A       ENSG00000148303  2          4           ENST00000291527_0:551_ENST00000463740_29:1164;ENST00000291527_0:551_ENST00000323345_33:891
ACTG1       ENSG00000184009  CANX        ENSG00000127022  2          2           ENST00000575994_0:575_ENST00000505090_429:485;ENST00000575659_0:463_ENST00000505090_429:485;...
FAM177B     ENSG00000197520  SOD2        ENSG00000112096  2          4           ENST00000391880_0:251_ENST00000367054_632:815;ENST00000445590_0:251_ENST00000367054_632:815;...
EEF2        ENSG00000167658  S100A6      ENSG00000197956  2          2           ENST00000309311_0:3164_ENST00000462951_0:311;ENST00000600794_0:352_ENST00000462951_0:311
```

Converting pizzly. Exploring [pVAC code](https://github.com/griffithlab/pVACtools/blob/master/lib/pipeline.py#L115-L132). [`IntegrateConverter`](https://github.com/griffithlab/pVACtools/blob/935dbb31a92377dba5bc580db962dbfecc2be140/lib/input_file_converter.py#L333) and [`FusionFastaGenerator`](https://github.com/griffithlab/pVACtools/blob/935dbb31a92377dba5bc580db962dbfecc2be140/lib/fasta_generator.py#L212) are fusion specific, but only the former reads the INTEGRATE results; FastaGenerator seem to be working from internal intermediate files, so we should'n't care about it.

`IntegrateConverter` reads the following fields:

```
chromosome_name            : 'chr 5p' / 'chr 3p'
start                      : 'start 5p' / 'start 3p'
stop                       : 'end 5p' / 'end 3p'
gene_name                  : 'name of fusion'
protein_positions          : 'fusion positions'
transcript names           : 'transcripts' 
fusion_amino_acid_sequence : 'peptides'
variant_type               : 'can be in-frame' -> 'inframe_fusion' or 'frameshift_fusion'
```

Expected output (a header added based on the [source code](https://github.com/griffithlab/pVACtools/blob/935dbb31a92377dba5bc580db962dbfecc2be140/lib/input_file_converter.py#L336-L353)) is the following:

```
chr 5p  start 5p  end 5p    chr 3p  start 3p  end 3p    name of fusion  tier              can be in-frame  peptides                                                                          fusion positions     transcripts                                                                                                                
19      -1        39123136  6       -1        46607405  EIF3K>>CYP39A1  1     +  -  2  1  1                QALDENMDLLEGITGFEDSVRKSSIPKNVFLALHEKLYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLV  22                1  ENST00000248342(421)|ENST00000538434(160)|ENST00000545173(421)|ENST00000592558(421)|ENST00000593149(160);ENST00000275016;  1
17      -1        7490589   17      9760739   -1        MPDU1>>GLP2R    3     +  +  4  1  1                LLQAATNYHNGHTGQLSAITVFLLFGGSLARIFTSIQKTPLHAQLHPHELVCFFHPENPGCTGEGRRLLQLLLQEAX,SPGSHQLPQRAHRPALSHHSLPAVWGLPGPNLHFHSENSTARATTSTX,FSRQPPTTTTGTQASSQPSQSSCCLGAPWPESSLPFRKLHCTRNYIHMNLFASFILRTLAVLVKDVVFYNSYSKRPDNENGWMSYLSE  37,36,36  0,1,2  ENST00000250124(618);;ENST00000262441|ENST00000458005|ENST00000574745,ENST00000423172(532);;ENST00000262441|ENST00000458005|ENST00000574745,ENST00000396501(560);ENST00000262441|ENST00000574745;ENST00000458005  1
```

We can match the following fields:

```
chromosome_name            : -
start                      : -
stop                       : -
gene_name                  : 'geneA.name' >> 'geneB.name'
fusion_positions           : ?
transcript names           : 'transcripts.list' 
fusion_amino_acid_sequence : -
variant_type               : -
```

We don't have genomic coordinates, fusion positions, AA sequences, and if it's inframe.

### Oncofuse

Oncofuse runs only for GRCh37 projects because it needs STAR fusion calls from STAR alignment which doesn't work with hg38's alts.

On the bright side, Oncofuse outputs genomic coordinates:

```
SAMPLE_ID  FUSION_ID  TISSUE  SPANNING_READS  ENCOMPASSING_READS  GENOMIC                         5_FPG_GENE_NAME  5_IN_CDS?  5_SEGMENT_TYPE  5_SEGMENT_ID  5_COORD_IN_SEGMENT  5_FULL_AA  5_FRAME  3_FPG_GENE_NAME  3_IN_CDS?  3_SEGMENT_TYPE  3_SEGMENT_ID  3_COORD_IN_SEGMENT  3_FULL_AA  3_FRAME  FPG_FRAME_DIFFERENCE  P_VAL_CORR   DRIVER_PROB            EXPRESSION_GAIN       5_DOMAINS_RETAINED  3_DOMAINS_RETAINED \
           32         AVG     4               0                   chr11:125765392>chr11:3381400   PUS3             Yes        Exon            3             292                 223        0        ZNF195           Yes        Exon            6             395                 352        1        1                     1.0          0.6976964618223671     0.04827719710198153   Pseudouridine s
           90         AVG     4               0                   chr12:49957249>chr14:31085512   MCRS1            Yes        Exon            6             80                  225        0        G2E3             Yes        Exon            15            29                  77         1        1                     1.0          0.003852211235299983   NaN
           103        AVG     4               0                   chr21:43771051>chr6:74228281    TFF2             Yes        Exon            1             157                 4          2        EEF1A1           Yes        Exon            6             52                  189        2        1                     1.0          6.77529842030136E-5    -0.6484925527690931
  5_DOMAINS_BROKEN  3_DOMAINS_BROKEN  5_PII_RETAINED      3_PII_RETAINED  CTF                   G                     H                      K    P                      TF
```

Runs from STAR-Fusion output, which can be visualized in a circos with [fsnviz](https://github.com/bow/fsnviz).


### JAFFA

[https://github.com/Oshlack/JAFFA](https://github.com/Oshlack/JAFFA)

Does no output genome cooridnates either.


### Clinker

[Clinker](https://github.com/Oshlack/Clinker/wiki) can be used to map to genome reference and to visualize, and can work from Pizzly too. More details at the [publication](https://academic.oup.com/gigascience/advance-article/doi/10.1093/gigascience/giy079/5049009)


### SQUID

[https://github.com/Kingsford-Group/squid](https://github.com/Kingsford-Group/squid)

Can predict both fusions and non-fusion SVs from RNAseq.

Runs either from 2 STAR BAM files: alignment BAM file and chimeric BAM file:

```
squid -b alignment.bam -c chimeric.bam -o squidout
```

Or a combined BAM file of both concordant and discordant alignments generated by [BWA](http://bio-bwa.sourceforge.net/) or [SpeedSeq](https://github.com/hall-lab/speedseq):

```
squid --bwa -b combined_alignment.bam -o squidout
```

### INTEGRATE-Neo

We can try to annotate Pizzly calls with [INTEGRATE-NEO](https://github.com/ChrisMaherLab/INTEGRATE-Neo), and feed the result into pVACfuse.

Installing (spartan)

```
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/INTEGRATE-NEO
git clone https://github.com/ChrisMaherLab/INTEGRATE-Neo
cd INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1
module load CMake  # or `conda install cmake`
module load GCC/6.4.0-2.28
./install.sh -o /data/cephfs/punim0010/extras/vlad/miniconda/envs/fusions/bin
```

Running.

`integrate-neo.py` script is the whole pipeline and need fastq on input to call HLA. We don't need that, so will try to run it step by step.

[integrate-neo.jpg](integrate-neo.jpg)

#### fusionBedpeAnnotator

Input is technically a bedpe file produced by INTEGRATE, by can simulate it from pizzly. Needs to satisfy standards of The ICGC-TCGA DREAM Somatic Mutation Calling-RNA Challenge (SMC-RNA)- basically it's just a minimal bedpe file.

```
fusionBedpeAnnotator \
    --reference-file /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
    --gene-annotation-file /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file ./INTEGRATE-Neo/Examples/example1/fusions.bedpe \
    --output-file fusions.annotated.bedpe
```

Input:

```
12         -1   14456634  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -    1
12         -1   14457295  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -   42
21   41497936         -1  12   14460493         -1  TMPRSS2>>ATF7IP        3  -  +   16
21   41498118         -1  12   14460494         -1  TMPRSS2>>ATF7IP        3  -  +  214
```

Output:

```
2          -1   14456634  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -    1  1  0  EHFEYRCHFLCLIWLYVRRNLASFQVX  24  2   ENST00000261168(2069)|ENST00000536444(2066)|ENST00000540793(2069)|ENST00000543189(2066)|ENST00000544627(2093);;ENST00000332149|ENST00000398585|ENST00000458356|ENST00000463138                      0
12         -1   14457295  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -   42  1  1  NSKKVGKNRKCNRX               29  1   ENST00000261168(2158)|ENST00000536444(2155)|ENST00000540793(2158)|ENST00000543189(2155)|ENST00000544627(2182);ENST00000332149|ENST00000398585|ENST00000458356;ENST00000463138                       0
21   41497936         -1  12   14460493         -1  TMPRSS2>>ATF7IP        3  -  +   16  0  0  NA                           NA  NA  NA                                                                                                                                                                                                  0
21   41498118         -1  12   14460494         -1  TMPRSS2>>ATF7IP        3  -  +  214  1  0  X                            23  0   ENST00000332149(3055)|ENST00000398585(3053)|ENST00000458356(1597);;ENST00000261168|ENST00000536444|ENST00000537653|ENST00000539659|ENST00000540793|ENST00000541654|ENST00000543189|ENST00000544627  0
```

For some reason peptide sequences are different from the example output at `INTEGRATE-Neo/Examples/example1/fusion_antigen_out/fusions.bedpe.annot`, and the last column is all 0.

Maybe try to regenerate the `.genePred` file according to [pVAC docs](https://pvactools.readthedocs.io/en/latest/pvacfuse/prerequisites.html)?

```
# source ~/load_bcbio.sh
gtfToGenePred -genePredExt -geneNameAsName2 /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf hg38.genePred

fusionBedpeAnnotator \
    --reference-file /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
    --gene-annotation-file hg38.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file ./INTEGRATE-Neo/Examples/example1/fusions.bedpe \
    --output-file fusions.annotated.bedpe
```

Same results. Maybe it has to do with repeated error messages?

```
Inquire a not existed reference sequence in chrNameToIndex. exception caught: map::at
```

Trying to follow the [INTEGRATE-Neo docs](https://github.com/ChrisMaherLab/INTEGRATE-Neo) to create the .genePred file. Downloading the GTF file from Ensmebl and running gtfToGenePred on it:

```
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.86.gtf.gz Homo_sapiens.GRCh38.86.genePred

fusionBedpeAnnotator \
    --reference-file /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
    --gene-annotation-file Homo_sapiens.GRCh38.86.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file ./INTEGRATE-Neo/Examples/example1/fusions.bedpe \
    --output-file fusions.annotated.bedpe
```

Same result. Maybe try to download the Ensembl reference sequence as well?

```
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

fusionBedpeAnnotator \
    --reference-file Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --gene-annotation-file Homo_sapiens.GRCh38.86.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file ./INTEGRATE-Neo/Examples/example1/fusions.bedpe \
    --output-file fusions.annotated.bedpe
```

Works great now, the results are identical:

```
12         -1   14456634  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -    1  1  0  HPPNPPVSPGKTVNDVNSNNNMSYRDX                                 24  2   ENST00000261168(2069)|ENST00000536444(2066)|ENST00000540793(2069)|ENST00000543189(2066)|ENST00000544627(2093);;ENST00000332149|ENST00000398585|ENST00000458356|ENST00000463138                       1
12         -1   14457295  21         -1   41489593  ATF7IP>>TMPRSS2        3  +  -   42  1  1  NAGTVRQMLESKRNVSESAPPSFQTPVNTETKKALCITLTLGTFLVGAALAAGLLWKF  29  1   ENST00000261168(2158)|ENST00000536444(2155)|ENST00000540793(2158)|ENST00000543189(2155)|ENST00000544627(2182);ENST00000332149|ENST00000398585|ENST00000458356;ENST00000463138                        1
21   41497936         -1  12   14460493         -1  TMPRSS2>>ATF7IP        3  -  +   16  0  0  NA                                                          NA  NA  NA                                                                                                                                                                                                   0
21   41498118         -1  12   14460494         -1  TMPRSS2>>ATF7IP        3  -  +  214  1  0  HIEHSRYLSLLDAVDNSKMALNSYLQPILSLLQQLSVVNLNCRLQX              23  0   ENST00000332149(3055)|ENST00000398585(3053)|ENST00000458356(1597);;ENST00000261168|ENST00000536444|ENST00000537653|ENST00000539659|ENST00000540793|ENST00000541654|ENST00000543189|ENST00000544627   1
```

#### fusionBedpeSubsetter

Filter based on rules.

```
cat INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt
# rule    $(NF-4)!="NA" && $(NF-6)==1

fusionBedpeSubsetter --input-file fusions.annotated.bedpe --rule-file INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt --output-file fusions.annotated.filt.bedpe
```

Getting down to 18 from 26 rows.

#### runNetMHC4WithSMCRNABedpe.py

Predicts affinity. Needs HLA alleles on input, taking example HLAminer output at `./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv`

Peptide length: pVACfuse uses 21 as default, but INTEGRATE-Neo uses 8,9,10,11. Picking 8,9,10,11 for now.

Avalilable alleles: not clear what's this file, but trying pVACfuse's valid_alleles:

```
pvacfuse valid_alleles > pvacfuse_alleles.txt

runNetMHC4WithSMCRNABedpe.py \
    -o MHC \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f fusions.annotated.filt.bedpe \
    -p 8,9,10,11 \
    -v pvacfuse_alleles.txt \
    -n ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/vendor/netMHC-4.0/netMHC \
    -k
```

Fails with (in `MHC/netMHC4.0.out.append.txt`):

```
netMHC: no binaries found for Linux_x86_64 /usr/cbs/packages/netMHC/4.0/netMHC-4.0/Linux_x86_64/bin/netMHC
```

Copying MHC from pVAC isntallation from raijin:

```
scp -r rjn:/g/data3/gx8/projects/Saveliev_pVACtools/mhc_i/method/netmhc-4.0-executable/netmhc_4_0_executable .
```

Trying with it:

```
runNetMHC4WithSMCRNABedpe.py \
    -o MHC \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f fusions.annotated.filt.bedpe \
    -p 8,9,10,11 \
    -v pvacfuse_alleles.txt \
    -n ./netmhc_4_0_executable/netMHC \
    -k
```

Having issues with HLA against valid list. In fact we should use the list from MHC installation: `/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/INTEGRATE-NEO/netmhc_4_0_executable/Linux_x86_64/data/allelelist`

Also doing some fixes in the code ([INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/runNetMHC4WithSMCRNABedpe.py](INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/runNetMHC4WithSMCRNABedpe.py))

```
runNetMHC4WithSMCRNABedpe.py \
    -o MHC \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f fusions.annotated.filt.bedpe \
    -p 8,9,10,11 \
    -v ./netmhc_4_0_executable/Linux_x86_64/data/allelelist \
    -n ./netmhc_4_0_executable/netMHC \
    -k
```

Runs well, but interestingly only one HLA allele match the list (HLA-C*07:02).

#### runAddNetMHC4Result.py

```
runAddNetMHC4Result.py \
    -b fusions.annotated.filt.bedpe \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f ./MHC/netMHC4.0.out.append.txt \
    -o fusions.annotated.filt.mhc.bedpe
```

The result is 2 calls:

```
21      41498118        -1      12      14460494        -1      TMPRSS2>>ATF7IP 3       -       +       214     SYLQPILSL       256.56  HLA-C07:02      Top_1   4576    4.34e-79        783.6
1       225999192       -1      1       116019581       -1      SDE2>>SLC22A15  3       -       +       1       HRHCQDQWF       141.96  HLA-C07:02      Top_1   4576    4.34e-79        783.6
```

Only the 2nd one is reproted in the example data `./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/result.bedpe`, might be due to the default affinity score threshold of 500 (set with `-m` in `runAddNetMHC4Result.py`). If we run with -m 200, the result will be identical.


## Running pVACfuse with INTEGRATE-Neo output

(`raijin`)

Taking the HLA alleles from the example output `./INTEGRATE-NEO/INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv`: "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61"

```
cd /g/data3/gx8/projects/Saveliev_pVACtools/fusions
scp -r spa:/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/INTEGRATE-NEO/fusions.annotated.filt.mhc.bedpe INTEGRATE-NEO

pvacfuse run \
    INTEGRATE-NEO/fusions.annotated.filt.mhc.bedpe \
    Test \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NetMHC \
    pvacfuse_integrate_neo \
    -e 8,9,10,11 \
    --top-score-metric=lowest \
    --keep-tmp-files

(five_p_transcripts, three_p_inframe_transcripts, three_p_frameshift_transcripts) = transcript_set.split(';')
ValueError: not enough values to unpack (expected 3, got 1)
```

The error is in parsing the transcript column, as it's missing transcripts from the expected format:

```
chr 5p  start 5p  end 5p    chr 3p  start 3p  end 3p    name of fusion  tier              can be in-frame  peptides                                                                          fusion positions     transcripts                                                                                                                
19      -1        39123136  6       -1        46607405  EIF3K>>CYP39A1  1     +  -  2  1  1                QALDENMDLLEGITGFEDSVRKSSIPKNVFLALHEKLYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLV  22                1  ENST00000248342(421)|ENST00000538434(160)|ENST00000545173(421)|ENST00000592558(421)|ENST00000593149(160);ENST00000275016;  1
```

The annotated `.bedpe` file (pre-MHC) had transcript annotations, so it might be a bug quick to fix:

```
scp -r spa:/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/INTEGRATE-NEO/fusions.annotated.bedpe INTEGRATE-NEO
grep 225999192 INTEGRATE-NEO/*
INTEGRATE-NEO/fusions.annotated.bedpe:1           225999192  -1  1  116019581  -1  SDE2>>SLC22A15  3  -  +  1  1          1       MAEAAALVWIRGPGFGCKAVRCASGRCTVRDFIHRHCQDQWFLIANRSYKVSAASSFFFSGVFVGVISFGQLSDRFGRKKVYLT  40        0  ENST00000272091(3836);ENST00000369502|ENST00000369503;  1
INTEGRATE-NEO/fusions.annotated.filt.mhc.bedpe:1  225999192  -1  1  116019581  -1  SDE2>>SLC22A15  3  -  +  1  HRHCQDQWF  141.96  HLA-C07:02                                                                            Top_1  4576  4.34e-79                                                783.6
```

We exprect the following file:

```
chr 5p  start 5p   end 5p  chr 3p  start 3p   end 3p  name of fusion  tier              can be in-frame  peptides                                                                              fusion positions     transcripts                                                                                                                
1       225999192  -1      1       116019581  -1      SDE2>>SLC22A15  3     -  +  1  1  1                MAEAAALVWIRGPGFGCKAVRCASGRCTVRDFIHRHCQDQWFLIANRSYKVSAASSFFFSGVFVGVISFGQLSDRFGRKKVYLT  40                0  ENST00000272091(3836);ENST00000369502|ENST00000369503;  1
```

Which matches perfectly the annotated, but pre-netMHC. That makes sence as pVACfuse re-runs netMHC. So using the annotated file (or even subsetted?):

```
scp -r spa:/data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/INTEGRATE-NEO/fusions.annotated.filt.bedpe INTEGRATE-NEO

pvacfuse run \
    INTEGRATE-NEO/fusions.annotated.filt.bedpe \
    Test \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NetMHC \
    pvacfuse_integrate_neo \
    -e 8,9,10,11 \
    --top-score-metric=lowest \
    --keep-tmp-files
```

Works great, though again  Outputs 2 epitopes:

```
Chromosome  Start                  Stop     Reference  Variant  Ensembl Gene ID  Variant Type    Mutation  Protein Position  Gene Name        HLA Allele   Peptide Length  Sub-peptide Position  Mutation Position  MT Epitope Seq  WT Epitope Seq  Best MT Score Method  Best MT Score  Corresponding WT Score  Corresponding Fold Change  Tumor DNA Depth  Tumor DNA VAF  Tumor RNA Depth  Tumor RNA VAF  Normal Depth  Normal VAF  Gene Expression  Transcript Expression  Median MT Score  Median WT Score  Median Fold Change  NetMHC WT Score  NetMHC MT Scor
1 / 1       225999192 / 116019581  -1 / -1  fusion     fusion   NA               inframe_fusion  NA        40                SDE2>>SLC22A15   HLA-C*07:02  9               4                     NA                 HRHCQDQWF       NA              NetMHC                141.96         NA                      NA                         NA               NA             NA               NA             NA            NA          NA               NA                     141.96           NA               NA                  NA               141.96
21 / 12     41498118 / 14460494    -1 / -1  fusion     fusion   NA               inframe_fusion  NA        23                TMPRSS2>>ATF7IP  HLA-C*07:02  9               10                    NA                 SYLQPILSL       NA              NetMHC                256.56         NA                      NA                         NA               NA             NA               NA             NA            NA          NA               NA                     256.56           NA               NA                  NA               256.56
```

And they are the same as produced INTEGRATE-Neo itself.


## Converting pizzly into INTEGRATE-type bedpe

We need to convert pizzly output into a bedpe file. Basically we need to query coordinates for the Ensembl transcript IDs, and using them (plus strand sign) convert the transcript-based fusion coordinates into the genome-based coordinate.

Pizzly output looks as follows. Flat:

```
geneA.name  geneA.id         geneB.name  geneB.id         paircount  splitcount  transcripts.list
TFF1        ENSG00000160182  RPL7A       ENSG00000148303  2          4           ENST00000291527_0:551_ENST00000463740_29:1164;ENST00000291527_0:551_ENST00000323345_33:891
```

Corresponding json record:

```
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
...
```

We expect it to convert into following:

```
#chr 5p  start 5p  end 5p    chr 3p  start 3p   end 3p    name of fusion    tier  strand 3p  strand 5p  quantitation
21       -1        43786519  9       -1         136215101 TFF1>>RPL7A	    -     -          +          -
```

Making a script that maps transcript coordinates in pizzly output using [pyensembl](https://github.com/openvax/pyensembl) and produces a bedpe file.

```
# takes CCR180109_VPT-EH09_RNA.json and CCR180109_VPT-EH09_RNA-flat-filtered.tsv
# produces CCR180109_VPT-EH09_RNA.bedpe
python pizzly_to_bedpe.py CCR180109_VPT-EH09_RNA
```

We get:

```
21      -1        43786519        9       136215101       -1       TFF1>>RPL7A     3       -       +       6
...
```

Running INTEGRATE-Neo to annotate the bedpe. Need to make sure we use the same Ensembl version (75 for GRCh37 in this case).

```
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh37.75.gtf.gz Homo_sapiens.GRCh37.75.genePred

wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz

fusionBedpeAnnotator \
    --reference-file Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --gene-annotation-file Homo_sapiens.GRCh37.75.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file /data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/pizzly/CCR180081_MH18T002P053_RNA-flat-filtered.bedpe \
    --output-file /data/cephfs/punim0010/projects/Saveliev_pVACtools/fusion/pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe
```

All good. However, we are getting a lot of repetitive records when same fusion hits different transcripts. It's irrelevant for epitope prediction when all transcript fusions produce the same peptides, so should do some filtering based on that.

## pVACfuse on annotated pizzly 

Now running `pVACfuse` (`raijin`):

```
pvacfuse run \
    CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe \
    CCR180081_MH18T002P053_RNA \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NetMHC \
    pvacfuse_integrate_neo_pizzly \
    -e 8,9,10,11 \
    --top-score-metric=lowest \
    --keep-tmp-files
```

Doesn't output anything. On the other hand, only one allele HLA-C*07:02 is supported by NetMHC. Trying to add other tools, identical to the pVACseq:

```
pvacfuse run \
    CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe \
    CCR180081_MH18T002P053_RNA \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
    pvacfuse_integrate_neo_pizzly_all_callers \
    -e 8,9,10,11 \
    --top-score-metric=lowest \
    --keep-tmp-files
```

This outputs 337 epitopes in 9 gene fusions, 11 transcript fusions. All are in-frame. The majority is 9-peptide length, however a few are 8. 


## INTEGRATE-Neo on pizzly and adding filtering before pVACfuse:

Also running the remaining INTEGRATE-Neo pipeline to compare its results with NetMHC (`spartan`):
 
```
cat INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt
# rule    $(NF-4)!="NA" && $(NF-6)==1
fusionBedpeSubsetter \
    --input-file ../pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe \
    --rule-file INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt \
    --output-file ../pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe
# 142(45 unique) -> 53(11 unique)

runNetMHC4WithSMCRNABedpe.py \
    -o MHC \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f ../pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe \
    -p 8,9,10,11 \
    -v ./netmhc_4_0_executable/Linux_x86_64/data/allelelist \
    -n ./netmhc_4_0_executable/netMHC \
    -k

runAddNetMHC4Result.py \
    -b ../pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe \
    -a ./INTEGRATE-Neo/Examples/example1/fusion_antigen_out/HLAminer_alleles.tsv \
    -f ./MHC/netMHC4.0.out.append.txt \
    -o ../pizzly/CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt-mhc.bedpe
```

Empty output. Which aligns with an empty output from pVACseq when running NetMHC only.

However, we can steal the subsetting method to perform before pVACfuse run (`raijin`):

```
pvacfuse run \
    CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe \
    CCR180081_MH18T002P053_RNA \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
    pvacfuse_integrate_neo_pizzly_all_callers_subset \
    -e 8,9 \
    --top-score-metric=lowest \
    --keep-tmp-files
```

Getting same 337 epitopes, same as before filtering.









