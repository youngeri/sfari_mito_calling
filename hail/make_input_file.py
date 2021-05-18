#!/usr/bin/env python

'''
make input files that contain lists of paths, samples
check if sany samples have '.' in name!
'''


# test run of 112 variant calling  control samples
test_dir = '/media/erica/Seagate/research/callset_20200225/mt-data/test/'
vcf_only_file = test_dir + 'vcf_only_list.txt'
coverage_only_file = test_dir + 'coverage_only_list.txt'

# to write: tab delimited file with sample, base_level_coverage_metrics path
coverage_list_file = test_dir + 'input_coverage_file.tsv'

# to write: tab delimited file with sample, vcf path
vcf_list_file = test_dir + 'input_vcf_file.tsv'

sample_dict = {}


with open(vcf_list_file, 'w') as vcf_list:
	with open(vcf_only_file, 'r') as vcf_only:
		for vcf_line in vcf_only:
			vcf_path = vcf_line.strip()
			sample_id = vcf_line.split('/')[-1].replace('.final.split.vcf', '').strip()
			
			run_number = vcf_line.split('/')[6]
			
			sample_dict[run_number] = sample_id
			vcf_list.write('\t'.join([sample_id, vcf_path]) + '\n')


with open(coverage_list_file, 'w') as coverage_list:
	with open(coverage_only_file, 'r') as coverage_only:
		for coverage_line in coverage_only:
			coverage_path = coverage_line.strip()
			run_number = coverage_path.split('/')[6]
			sample_id = sample_dict[run_number]
			coverage_list.write('\t'.join([sample_id, coverage_path]) + '\n')
