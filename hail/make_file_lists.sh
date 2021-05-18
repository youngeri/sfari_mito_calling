#!/bin/bash

# make lists of VCF and coverage files
# some files get stored in attempt2 directories, some do not

GCP_DIR="gs://wustl-mt-pipeline-202010/vcc/outputs/MitochondriaPipeline" 
LOCAL_DIR="/media/erica/Seagate/research/callset_20200225/mt-data/test"

gsutil ls -r $GCP_DIR | grep final.split.vcf > $LOCAL_DIR/vcf_only_list.txt
gsutil ls -r $GCP_DIR | grep per_base_coverage.tsv > $LOCAL_DIR/coverage_only_list.txt
 
