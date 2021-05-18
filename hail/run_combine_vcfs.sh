#!/bin/bash

# test run
DIR="gs://epytest/data/mt-202010/test"

hailctl dataproc submit epy-mt \
combine_vcfs.py \
--coverage_mt_path $DIR/merged_coverage.mt \
--input_tsv $DIR/input_vcf_file.tsv \
--output_bucket $DIR \
--file_suffix "test_sequencing_controls" \
