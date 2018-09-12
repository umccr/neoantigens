# KRAS-wt patient (MON3)

- [KRAS-wt patient (MON3)](#kras-wt-patient-mon3)
  - [Data (raijin):](#data-raijin)
  - [Preparing data](#preparing-data)
  - [Running pVACseq](#running-pvacseq)
  - [Expression](#expression)
- [TODO](#todo)

This patient has over 200k somatic mutations.

## Data (raijin):

```
# Tumor BAM
/g/data/gx8/data/pVAC/hg38_wgs_samples/final/PRJ170194_GBB10_T/PRJ170194_GBB10_T-ready.bam

# VCF
/g/data/gx8/data/pVAC/hg38_wgs_samples/final/2018-09-11_2018-09-04T0624_pvac_exploration_WGS-merged/MON3-ensemble-annotated.vcf.gz

# Tumor HLA
/g/data3/gx8/data/pVAC/hg38_wgs_samples/final/PRJ170194_GBB10_T/PRJ170194_GBB10_T-hla-optitype.csv
# Normal HLA
/g/data3/gx8/data/pVAC/hg38_wgs_samples/final/PRJ170195_GBB10_B/PRJ170195_GBB10_B-hla-optitype.csv
```

## Preparing data

(`Raijin`)

Working at `/g/data3/gx8/projects/Saveliev_pVACtools/MON3`

```
# VCF: filtering with panel of normals
pon_anno data/MON3-ensemble-annotated.vcf.gz -h1 | bcftools view -f.,PASS -o data/MON3-ensemble-pon_anno.h1.pass.vcf

# VCF: annotating
vep \
--input_file data/MON3-ensemble-pon_anno.h1.pass.vcf --format vcf \
--output_file data/MON3-ensemble-pon_anno.h1.pass.VEP.vcf --vcf \
--pick --symbol --terms SO --plugin Downstream --plugin Wildtype \
--cache --dir_cache ../vep_data/GRCh38 --assembly GRCh38 --offline

# VCF: getting bcftools and extracting the tumor sample:
bcftools view -s PRJ180359_MON3-T data/MON3-ensemble-pon_anno.h1.pass.VEP.vcf > data/MON3-ensemble-pon_anno.h1.pass.VEP.TUMOR.vcf

# HLA - same for tumor and normal:
HLA-A*03:01,HLA-A*26:01,HLA-B*55:01,HLA-B*38:01,HLA-C*03:03,HLA-C*12:03
```

## Running pVACseq

```
VCF=data/MON3-ensemble-pon_anno.h1.pass.VEP.TUMOR.vcf
SAMPLE=PRJ180359_MON3-T
HLA_TYPES="HLA-A*03:01,HLA-A*26:01,HLA-B*55:01,HLA-B*38:01,HLA-C*03:03,HLA-C*12:03"

pvacseq run \
$VCF \
$SAMPLE \
"$HLA_TYPES" \
NetMHCIIpan NetMHCcons SMM SMMPMBEC NNalign SMMalign \
pvacseq_results \
-e 8,9,10,11 \
-t \
--top-score-metric=lowest \
--iedb-install-directory /g/data3/gx8/projects/Saveliev_pVACtools \
--trna-vaf 10 \
--net-chop-method cterm \
--netmhc-stab \
--exclude-NAs
```

## Expression
# TODO


Since we're doing all sensible somatic filtering before passing the VCF into pVACtools, we don't need to do extra filtering here. For example, filtering out normal AF>2% will have trouble with tumor contamination, and want to deal with that before running pVAC. So we will just remove the frequencies and depths from the inputs (except for those for RNA).

However, keeping in mind that we want to use only mutations for building epitopes that consistently appear in a tumor altogether. It's reached ideally on the VAF roughly equal to tumor ploidy. That's why it's good to figure out tumor ploidiy based on typical somatic mutations AFs, and then filter out mutations with tumor VAF below that value.

Running coverage to consider it later:

(`Spartan`)

Splitting the variants into SNPs and indels, and converting to site-list files for bam-readcount:

```
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt
mkdir coverage
VCF=MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf
bcftools view -v snps $VCF | bcftools query -f "%CHROM\t%POS\n" | awk '{OFS="\t"}{ print $1,$2,$2 }' > coverage/tumor_snps_sites.txt
bcftools view -v indels $VCF | bcftools query -f "%CHROM \t %POS \t %LEN \n" | awk '{OFS="\t"}{ print $1,$2,$2+$3-1 }' > coverage/tumor_indels_sites.txt
```

Installing bam-readcount with conda into the `pvac` environment causes a conflict with `perl-storable`, thus creating a new environment with `conda create -n bam-readcount bam-readcount`.

Running bam-readcount for the NeverResponder (`-w 1` to suppress repeated warnings, `-i` for indels).

```
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt
T_BAM=PRJ180359_MON3-T-ready.bam
N_BAM=PRJ180351_MON3-B-ready.bam
FA=/data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa

bam-readcount -f $FA -w 1 -l coverage/tumor_snps_sites.txt $T_BAM > coverage/tumor_snps.readcount &
bam-readcount -f $FA -w 1 -l coverage/tumor_indels_sites.txt -i $T_BAM > coverage/tumor_indels.readcount &
bam-readcount -f $FA -w 1 -l coverage/tumor_snps_sites.txt $N_BAM -i > coverage/normal_snps.readcount &
bam-readcount -f $FA -w 1 -l coverage/tumor_indels_sites.txt $N_BAM -i > coverage/normal_indels.readcount &
```

Creating yaml file `additional_input_file_list.yaml`:

```
normal_snvs_coverage_file:   /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt/coverage/normal_snps.readcount
normal_indels_coverage_file: /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt/coverage/normal_indels.readcount
tdna_snvs_coverage_file:     /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt/coverage/tumor_snps.readcount
tdna_indels_coverage_file:   /data/cephfs/punim0010/projects/Saveliev_pVACtools/KRAS-wt/coverage/tumor_indels.readcount
```