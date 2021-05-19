#!/bin/bash

HAPLOGREP="/home/erica/software/haplogrep/haplogrep"
DIR="/media/erica/Seagate/research/callset_20200225/mt-data/test"
INPUT=$DIR/combined_test_sequencing_controls_withGT.vcf.bgz
OUTPUT=$DIR/combined_test_sequencing_controls_haplogroups.txt

echo $INPUT

$HAPLOGREP classify \
--in $INPUT \
--format vcf \
--out $OUTPUT
