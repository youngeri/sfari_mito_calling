#!/bin/bash

# test run
DIR="gs://epytest/data/mt-202010/test"

hailctl dataproc submit epy-mt \
annotate_coverage.py \
--input_tsv ${DIR}/input_coverage_file.tsv \
--output_ht ${DIR}/merged_coverage.ht
