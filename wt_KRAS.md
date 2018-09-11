# KRAS-wt patient

- [KRAS-wt patient](#kras-wt-patient)
  - [VCF:](#vcf)
  - [HLA](#hla)
  - [Running pVACseq](#running-pvacseq)
  - [Coverage](#coverage)

This patient has over 200k somatic mutations.

Data on spartan:

```
/data/cephfs/punim0010/data/Results/Research-Croagh-KRAS-wt/2018-08-18/umccrised/MON3__PRJ180359_MON3-T/small_variants/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.vcf.gz
```

Copying data to Raijin into `/g/data3/gx8/projects/Saveliev_pVACtools/data/kras_wt` 

## VCF:

(`Raijin`)

```
vep \
--input_file data/kras_wt/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.vcf.gz \
--format vcf \
--pick \
--output_file data/kras_wt/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.VEP.vcf \
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
bcftools view -s PRJ180359_MON3-T data/kras_wt/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.VEP.vcf > data/kras_wt/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf
```

## HLA

(`Spartan`)

Realinging to hg38 HLA chromosomes. Build HLA-only reference:

```
cd /data/cephfs/punim0010/projects/Saveliev_pVACtools/hla_reference
samtools faidx /data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa $(cut -f1 /data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa.fai | grep HLA) > hg38_hla.fa
bwa index hg38_hla.fa
```

Getting bwa-kit to re-align:

```
conda install -c bioconda -y bwakit samtools htsbox  # samtools and htsbox executables must be in the same path as bwakit executable
```

Running (-S to disable BAM shuffling, -H to make HLA typing, -k to keep tmp files from HLA typing, -o is prefix)

```
mkdir hla_typing
run-bwamem -S -H -k -o hla_typing/PRJ180359_MON3-T_bwakit -t 30 ../hla_reference/hg38_hla.fa PRJ180359_MON3-T-ready.bam
# Will output a shell command. Replace `cat .bam` to avoid EOF errors:
htsbox bam2fq -O -t PRJ180359_MON3-T-ready.bam \
  | bwa mem -p -t30 -H'@RG\tID:PRJ180359_MON3-T\tSM:PRJ180359_MON3-T\tPL:illumina\tPU:PRJ180359_MON3-T' -C ../hla_reference/hg38_hla.fa - 2> hla_typing/PRJ180359_MON3-T_bwakit.log.bwamem \
  | samtools view -1 - > hla_typing/PRJ180359_MON3-T_bwakit.aln.bam
```

## Running pVACseq

```
VCF=data/kras_wt/MON3__PRJ180359_MON3-T-somatic-ensemble-pon_hardfiltered.VEP.TUMOR.vcf
SAMPLE=PRJ180359_MON3-T
HLA_TYPES="HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01"

pvacseq run \
$VCF \
$SAMPLE \
"$HLA_TYPES" \
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

## Coverage

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