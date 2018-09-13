# vaxrank

[VaxRank](https://github.com/openvax/vaxrank) is enother epitope prediction tool.

* Lack of docs (can refer to https://github.com/openvax/neoantigen-vaccine-pipeline/blob/master/README.md)
* Uses only one predictor at a time (though have a large choise: mhcflurry,netmhc,netmhc3,netmhc4,netmhccons,netmhccons-iedb,netmhciipan,netmhciipan-iedb,netmhcpan,netmhcpan-iedb,netmhcpan28,netmhcpan3,random,smm-iedb,smm-pmbec-iedb)
* Can't use fusions
* Runs from VCF and a BAM
* HLAs are expected in a weird notation

## Installation

```
pip install vaxrank
pyensembl install --release 75 --species human
```

## Running

First trying with netmhccons. Also can use NetMHCIIpan as it output something in pVACseq.

```
mkdir output_netmhccons
vaxrank \
    --vcf /data/cephfs/punim0010/data/Results/Patients/2018-01-17/umccrised/diploid__diploid_tumor/small_variants/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.vcf.gz \
    --genome GRCh37 \
    --bam /data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam \
    --vaccine-peptide-length 25 \
    --mhc-predictor netmhccons \
    --mhc-alleles "HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01" \
    --padding-around-mutation 5 \
    --output-ascii-report output_netmhccons/vaxrank-peptides.txt \
    --output-html-report output_netmhccons/vaxrank-peptides.html \
    --output-neoepitope-report output_netmhccons/vaxrank-neoepitope.txt \
    --log-path output_netmhccons/vaxrank.log
```

Doesn't run: `FileNotFoundError: [Errno 2] No such file or directory: 'netMHCcons': 'netMHCcons'`. 

Trying with `netmhccons-iedb` so it uses the IEDB instead of running the tool? 

No, needs internet to access the iedb database. No option to point to a local IEDB like pVACseq. Running on the login node.

Fails with exceptions

```
2018-09-12 19:23:40,046 - vaxrank.epitope_prediction:286 - INFO - 636 total peptides: 408 occur in reference, 578 failed score threshold
2018-09-12 19:23:40,063 - vaxrank.core_logic:191 - INFO - Keeping 8/11 vaccine peptides for Variant(contig='MT', start=8950, ref='G', alt='A', reference_name='GRCh37')
2018-09-12 19:36:07,146 - vaxrank.cli:326 - INFO - About to save args: {'vcf': ['/data/cephfs/punim0010/data/Results/Patients/2018-01-17/umccrised/diploid__diploid_tumor/small_variants/diploid__diploid_tumor-somatic-ensemble-pon_hardfiltered.vcf.gz'], 'maf': [], 'variant': [], 'genome': 'GRCh37', 'json_variants': [], 'bam': '/data/cephfs/punim0010/data/Results/Patients/2018-05-17/final/Unknown_B_RNA/Unknown_B_RNA-ready.bam', 'min_mapping_quality': 1, 'use_duplicate_reads': False, 'drop_secondary_alignments': False, 'min_alt_rna_reads': 2, 'min_variant_sequence_coverage': 2, 'variant_sequence_assembly': True, 'protein_sequence_length': 20, 'max_reference_transcript_mismatches': 2, 'include_mismatches_after_variant': False, 'min_transcript_prefix_length': 10, 'mhc_predictor': 'netmhccons-iedb', 'mhc_predictor_models_path': None, 'mhc_peptide_lengths': None, 'mhc_epitope_lengths': None, 'mhc_alleles_file': None, 'mhc_alleles': 'HLA-A*02:01,HLA-A*26:01,HLA-B*35:02,HLA-B*18:01,HLA-C*04:01,HLA-C*05:01', 'vaccine_peptide_length': 25, 'padding_around_mutation': 5, 'max_vaccine_peptides_per_mutation': 1, 'min_epitope_score': 1e-10, 'output_patient_id': '', 'output_csv': '', 'output_ascii_report': 'output_netmhccons/vaxrank-peptides.txt', 'output_html_report': 'output_netmhccons/vaxrank-peptides.html', 'output_pdf_report': '', 'output_json_file': '', 'output_xlsx_report': '', 'output_neoepitope_report': 'output_netmhccons/vaxrank-neoepitope.txt', 'num_epitopes_per_peptide': None, 'output_reviewed_by': '', 'output_final_review': '', 'log_path': 'output_netmhccons/vaxrank.log', 'max_mutations_in_report': None, 'output_passing_variants_csv': '', 'manufacturability': True, 'wt_epitopes': True, 'cosmic_vcf_filename': ''}
Traceback (most recent call last):
  File "/data/cephfs/punim0010/extras/vlad/miniconda/envs/fusions/bin/vaxrank", line 11, in <module>
    sys.exit(main())
  File "/data/cephfs/punim0010/extras/vlad/miniconda/envs/fusions/lib/python3.6/site-packages/vaxrank/cli.py", line 387, in main
    excel_report_path=args.output_neoepitope_report)
  File "/data/cephfs/punim0010/extras/vlad/miniconda/envs/fusions/lib/python3.6/site-packages/vaxrank/report.py", line 504, in make_minimal_neoepitope_report
    '%.2f nM' % epitope_prediction.wt_ic50),
TypeError: must be real number, not NoneType
```
