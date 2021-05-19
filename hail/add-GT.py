'''
add GT field to VCF for running haplogroup

if HL > 0.9, then GT = 1
else GT = 0
'''

import hail as hl

data_dir = 'gs://epytest/data/mt-202010/test/'
input_mt = data_dir + 'combined_test_sequencing_controls.mt'
output_vcf = data_dir + 'combined_test_sequencing_controls_withGT.vcf.bgz'

mt = hl.read_matrix_table(input_mt)
mt = mt.annotate_entries(GT = hl.if_else(mt.HL > 0.9, 1, 0))
hl.export_vcf(mt, output_vcf)
