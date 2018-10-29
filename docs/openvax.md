# OpenVax

## vaxrank

[VaxRank](https://github.com/openvax/vaxrank) is another wrapper around epitope predictors.

* Lack of docs. Can refer to https://github.com/openvax/neoantigen-vaccine-pipeline/blob/master/README.md, or to the [Snake rule](https://github.com/openvax/neoantigen-vaccine-pipeline/blob/master/pipeline/special_sauce.rules)
* Uses only one predictor at a time (though have a large choise: mhcflurry, netmhc, netmhc3, netmhc4, netmhccons(-iedb), netmhciipan(-iedb), netmhcpan(-iedb), netmhcpan28, netmhcpan3, random, smm-iedb, smm-pmbec-iedb)
* Can phase adjacent somatic/germline variants using a mutant coding sequence assembled from RNA reads
* Can't use fusions
* Runs from a VCF and a BAM
* HLAs are expected in a weird notation, however [mhcnames](https://github.com/openvax/mhcnames) can be helpful




### Installation

```
pip install vaxrank
pyensembl install --release 75 --species human
```


### Running

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


### Predictors

vaxrank uniquely supports mhcflurry, and pVACseq uniquely supports SMMalign and NNalign methods.

```
vaxrank                pVACtools                          MHC class     approach                                
------------------------------------------------------------------------------------------------------------------
netmhccons(,-iedb)     NetMHCcons                         I             consensus: NetMHC, NetMHCpan and PickPocket                         
netmhc(,3,4)           NetMHC     (part of NetMHCcons)    I             artificial neural networks (allele-specific)                                
netmhcpan(-iedb,3.28)  NetMHCpan  (part of NetMHCcons)    I             artificial neural networks (pan-specific)                                      
-                      PickPocket (part of NetMHCcons)    I             receptor-pocket similarity matrix, uses SMMPMBEC
netmhciipan(,-iedb)    NetMHCIIpan                        II            artificial neural networks (pan-specific)
smm-iedb               SMM                                I             similarity matrix
smm-pmbec-iedb         SMMPMBEC                           I             similarity matrix
-                      SMMalign                           II            stabilization matrix alignment       
-                      NNalign                            I, II         neural network models
mhcflurry              -                                  I             keras neural network               
```


## mhcflurry

[mhcflurry](https://github.com/openvax/mhcflurry) is a tool for preduction of peptides affinity in context of MHC I, that is not used in pVACtools, and is made by the OpenVax group (pyensemble, vaxrank).

### Installation (spartan; raijin doesn't work)

```
conda install -y numpy scipy pandas scikit-learn biopython
pip install mhctools mhcflurry
mhcflurry-downloads fetch
mhcflurry-downloads fetch data_iedb
```

### Running

On NeverResponder single fusion peptide:

```
mhcflurry-predict --alleles HLA-A0201 --peptides SLYKNKEV --out predictions.csv
```

Works!

```
tr ',' '\t' < predictions.csv | cols
allele     peptide   mhcflurry_prediction  mhcflurry_prediction_low  mhcflurry_prediction_high  mhcflurry_prediction_percentile
HLA-A0201  SLYKNKEV  1911.7888888097027    950.6672790559436         2852.644879789233          3.1505250000000014
```

Compared to pVACfuse:

```
Chromosome  Start    Stop                  Reference  Variant  Transcript                       Ensembl Gene ID  Variant Type    Mutation  Protein Position  Gene Name        HLA Allele   Peptide Length  Sub-peptide Position  Mutation Position  MT Epitope Seq  WT Epitope Seq  Best MT Score Method  Best MT Score 
12 / 6      -1 / -1  103932418 / 73520056  fusion     fusion   ENST00000614327-ENST00000309268  NA               inframe_fusion  NA        14                HSP90B1>>EEF1A1  HLA-A*02:01  8               1                     NA                 SLYKNKEV        NA              NetMHCcons            375.87        
12 / 6      -1 / -1  103932418 / 73520056  fusion     fusion   ENST00000614327-ENST00000309268  NA               inframe_fusion  NA        14                HSP90B1>>EEF1A1  HLA-A*02:01  8               2                     NA                 LYKNKEVS        NA              SMM                   419.962004535 
```


