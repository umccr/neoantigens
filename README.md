# 🏺 Neoantigen identification <!-- omit in toc --> 

Anti-tumor T cells recognize tumor somatic mutations, translated as single amino acid substitutions, as "neoantigens" (upon processing and presentation in the context of MHC molecules). From NGS somatic variant calling data, we can attempt to reconstruct peptide sequences of such neoantigens, and thus produce personalised cancer vaccines.

- [Introduction](#introduction)
    - [Biology](#biology)
    - [Approach](#approach)
    - [Prediction methods](#prediction-methods)
    - [Ranking tools](#ranking-tools)
- [pVACseq](#pvacseq)
- [wt_KRAS patient](#wtkras-patient)
- [pVACfuse](#pvacfuse)
- [NeoepitopePred](#neoepitopepred)
- [vaxrank](#vaxrank)

## Introduction

#### Biology

The main paradigm behind the development of cancer vaccines rests on the assumption that if the immune system is stimulated to recognize neoantigens, it may be possible to elicit the selective destruction of tumor cells. Vaccines incorporate these neoantigen peptides with the aim of enhancing the immune system’s anti-tumor activity by selectively increasing the frequency of specific CD8+ T cells, and hence expanding the immune system’s ability to recognize and destroy cancerous cells. This process is dependent on the ability of these peptides to bind and be presented by HLA class I molecules, a critical step to inducing an immune response and activating CD8+ T cells.

#### Approach

As we move from vaccines targeting ‘shared’ tumor antigens to a more ‘personalized’ medicine approach, in silico strategies are needed to first identify, then determine which somatic alterations provide the optimal neoantigens for the vaccine design. Ideally, an optimal strategy would intake mutation calls from tumor/normal DNA sequencing data, identify the neoantigens in the context of the patient’s HLA alleles, and parse out a list of optimal peptides for downstream testing.

#### Prediction methods

Several in silico epitope binding prediction methods have been developed [11](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=12175724), [12](http://scholar.google.com/scholar_lookup?title=A%20hybrid%20approach%20for%20predicting%20promiscuous%20MHC%20class%20I%20restricted%20T%20cell%20epitopes&author=M.%20Bhasin&author=G.%20Raghava&journal=J%20Biosci&volume=1&issue=32&pages=31-42&publication_year=2006), [13](https://scholar.google.com/scholar?hl=en&q=Lundegaard%20C%2C%20Lamberth%20K%2C%20Harndahl%20M%2C%20Buus%20S%2C%20Lund%20O%2C%20Nielsen%20M.%20NetMHC-3.0%3A%20accurate%20web%20accessible%20predictions%20of%20human%2C%20mouse%20and%20monkey%20MHC%20class%20I%20affinities%20for%20peptides%20of%20length%208-11.%20Nucleic%20Acids%20Res.%202008%3B36%28Web%20Server%20issue%29%3AW509%E2%80%93512.), [14](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=12717023), [15](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=25464113). These methods employ various computational approaches such as Artificial Neural Networks (ANN) and Support Vector Machines (SVM) and are trained on binding to different HLA class I alleles to effectively identify putative T cell epitopes.

#### Ranking tools

There are also existing tools ([IEDB](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=25300482), [EpiBot](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=25905908), [EpiToolKit](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=25712691)) that compile the results generated from individual epitope prediction algorithms to improve the prediction accuracy with consensus methods or a unified final ranking. EpiToolKit also has the added functionality of incorporating sequencing variants in its Galaxy-like epitope prediction workflow (via its Polymorphic Epitope Prediction plugin). However, it does not incorporate sequence read coverage or gene expression information available, nor can it compare the binding affinity of the peptide in the normal sample (WT) versus the tumor (mutant). Another multi-step workflow [Epi-Seq](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=25245761) uses only raw RNA-Seq tumor sample reads for variant calling and predicting tumor-specific expressed epitopes.

[pVACseq](docs/pVACseq.md) is a workflow that identifies and shortlist candidate neoantigen peptides from tumor mutational repertoire. It can use both DNA and RNA sequencing data. Predicted peptides can be used in a personalized vaccine after immunological screening. The tool offers the functionality to compare and differentiate the epitopes found in normal cells against the neoepitopes specifically present in tumor cells for use in personalized cancer vaccines.

[vaxrank](docs/vaxrank.md) is another workflow that prioritizes epitopes after running them agains prediction software like NetMHC.

[NeoepitopePred](#neoepitopepred) is an online-hosted workflow, can run from either single mutations or fusions.


## pVACseq

[pVACfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html) identifies and prioritizes neoantigens from a list of tumor mutations, and pre-calculated HLA types.

[Here we are exploring applitation of the tool](docs/pVACseq.md). In summary, we managed to run it from ensemble somatic variant calls generated by bcbio+umccrise, using HLA types predicted by Optitype in bcbio. The current pipeline is as follows.

```
# Run VEP to annotate the ensemble calls:
vep \
--input_file sample-somatic-ensemble-pon_hardfiltered.vcf.gz \
--format vcf \
--pick \
--output_file sample-somatic-ensemble-pon_hardfiltered.VEP.vcf \
--vcf --symbol --terms SO \
--plugin Downstream \
--plugin Wildtype \
--cache \
--dir_cache vep_data/GRCh37 \
--assembly GRCh37 \
--offline

# Subset the calls 
bcftools view -s sample sample-somatic-ensemble-pon_hardfiltered.VEP.vcf > sample-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf

# Take HLA predictions out of bcbio-nextgen Opitype output (hg38 only)

# Run pVACseq:
pvacseq run \
sample-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf \
sample \
"HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01" \
NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
pvacseq_results \
-e 9,10 \
-t \
--top-score-metric=lowest \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10 \
--net-chop-method cterm \
--netmhc-stab \
--exclude-NAs
```

For offline runs, we would omit `--net-chop-method cterm` and `--netmhc-stab`

Would keep `-t` to picking only the top epitope for a mutation, keeping in mind that we can go back to intermediate `diploid_tumor.filtered.coverage.tsv`, and refilter with `pvacseq top_score_filter`.

Also, with known tumor ploidy, it would be improtant to still filter the mutations with tumor VAF below the ploidy, to make sure the mutations for building epitopes appear in the tumor altogether. In that case, we would add the coverage files back into the inputs yamls and set the `--tdna-vaf` option.


## wt_KRAS patient

[Here we are exploring](docs/wt_KRAS.md) neoepitopes in wt_KRAS patient.


## pVACfuse

[pVACfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html) is a part of pVACtools package. It predicts neoantigens produced specifically by gene fusion events.

[Here we are exploring how can we run the tool in house](docs/pVACfuse.md). The tool works with input generated by a fusion caller `INTEGRATE`. Unfortuantely, it's not documented, the code base is not supported, and it doesn't install in our HPCs. Out of all other fusion callers, the most relevant for us is `pizzly`, as it shows a high performance, is integrated in bcbio-nextgen, and can work with both GRCh37 and hg38. The downside is that it runs from pseudo-alignments and does not provide genomic coordinates. We managed to convert `pizzly` calls to genomic coordinates, and annotate it to generate a `INTEGRATE`-like output that is get accepted by `pVACfuse`. 

The current pipeline is as below, given that we have a standard `pizzly` output `CCR180081_MH18T002P053_RNA{-flat-filtered.tsv,.json,.fusions.fasta}`.

```
# annotate and convert pizzly to bedpe format
python pizzly_to_bedpe.py CCR180081_MH18T002P053_RNA
# produces CCR180081_MH18T002P053_RNA-flat-filtered.bedpe

# get reference data for annotation (can sync with bcbio an `pyensembl` in the future)
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh37.75.gtf.gz Homo_sapiens.GRCh37.75.genePred
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz

# run `INTEGRATE-Neo`-bundled utility tool to annotate fusion junctions
fusionBedpeAnnotator \
    --reference-file Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --gene-annotation-file Homo_sapiens.GRCh37.75.genePred \
    --di-file ./INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/difile.txt \
    --input-file CCR180081_MH18T002P053_RNA-flat-filtered.bedpe \
    --output-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe

# remove low-quality
fusionBedpeSubsetter \
    --input-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated.bedpe \
    --rule-file INTEGRATE-Neo/INTEGRATE-Neo-V-1.2.1/src/rule.txt \
    --output-file CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe

# run pVACfuse
pvacfuse run \
    CCR180081_MH18T002P053_RNA-flat-filtered-annotated-filt.bedpe \
    CCR180081_MH18T002P053_RNA-VACfuse \
    "HLA-A*33:24,HLA-B*55:29,HLA-B*15:63,HLA-C*07:02,HLA-C*04:61" \
    NNalign NetMHCIIpan NetMHCcons SMM SMMPMBEC SMMalign \
    pvacfuse_out \
    -e 8,9 \
    --top-score-metric=lowest
```

## NeoepitopePred

Another tool that can run both from mutations and fusions is [NeoepitopePred](https://stjude.github.io/sjcloud-docs/guides/tools/neoepitope/).

Unfortunately, it requires to upload FASTQ or BAM data on the server, which might be a problem for us.

* Runs from fastq (calls optitype directly on fastq) or BAM (extacts HLA reads and unmapped reads, converts to fastq, then calls optitype)
* Can run EITHER from mutations, or fusions (like pVACtools)
* Need mutations in a specific format
* Uses NetMHCcons (consensus on a set of NetMHC methods).

## vaxrank

[VaxRank](https://github.com/openvax/vaxrank) is enother epitope prediction tool.

* Lack of docs (can refer to https://github.com/openvax/neoantigen-vaccine-pipeline/blob/master/README.md)
* Uses only one predictor at a time (though have a large choise: mhcflurry,netmhc,netmhc3,netmhc4,netmhccons,netmhccons-iedb,netmhciipan,netmhciipan-iedb,netmhcpan,netmhcpan-iedb,netmhcpan28,netmhcpan3,random,smm-iedb,smm-pmbec-iedb)
* Can't use fusions
* Runs from VCF and a BAM
* HLAs are expected in a weird notation

[Exploring it here](docs/vaxrank.md).