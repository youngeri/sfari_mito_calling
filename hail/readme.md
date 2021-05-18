# local data: /media/erica/Seagate/research/callset_20200225/mt-data/test
# GCP input data (gCVFs, coverage files): gs://wustl-mt-pipeline-202010/vcc/outputs/MitochondriaPipeline/
# GCP output data: gs://epytest/data/mt-202010/test

- make list of output VCFs and coverage files from mutect
  - issue is that some files get put in a directory called attempt2
  - script: make_file_lists.sh
  - output: vcf_only_list.txt, coverage_only_list.txt 

- make files that map sample ids and coverage files, vcfs
  - script: make_input_file.py
  - output: input_vcf_file.tsv, input_coverage_file.tsv
    - rsync to GCP

- make dataproc cluster (epy-mt)
  - script: create_hailctl

- combine coverage files
  - scripts: annotate_coverage.py, run_annotate_coverage.sh
  - input: input_coverage_file.tsv
  - output: 
    - merged_coverage.tsv, merged_coverage.ht - per site coverage statistics (mean, median, proportion over 100 or 1000)
    - merged_coverage_sample_level.txt - per sample coverage for each site
    - merged_coverage.mt - combination of data from above 2 files
    - temp directory

- merge VCFs
  - scripts: combine_vcfs.py, run_combine_vcfs.sh, raw_combined_mt.mt
  - output: combined_test_sequencing_controls.vcf.bgz, combined_test_sequencing_controls.mt

- make plots
  - gs://epytest/EOCAD_testscripts/sample_QC/mt-202010/trial-data-qc.ipynb
