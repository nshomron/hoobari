from stderr import printerr
import vcf
import pandas as pd
import re
import os
import subprocess
import vcfuid
import pprogress
import numpy as np

def str_to_int(gt):
	if gt is not None and gt is not '.':
		gt_split = gt.split('/')
		gt_sum = int(gt_split[0]) + int(gt_split[1])
		return gt_sum
	else:
		return None

def parse(vcf_file_path):
	printerr('parsing ' + vcf_file_path + '...')
	if os.path.isfile(vcf_file_path):
		index_length = pprogress.get_file_length(vcf_file_path)
		progress_index = pprogress.reset()

		vcf_list = []
		with open(vcf_file_path, 'r') as inp:
			reader = vcf.VCFReader(inp)

			samples = reader.samples
			vcf_list.append(['variant_name'] + samples)

			for record in reader:
				variant_name = [vcfuid.rec_to_uid(record)]
				genotypes = []
				for s in range(len(samples)):
					gt = str_to_int(record.genotype(samples[s]).data[0])
					genotypes.append(gt)
				line = variant_name + genotypes

				vcf_list.append(line)
				progress_index = pprogress.pprogress(progress_index, index_length)
			
			vcf_array = np.array(vcf_list)
			vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
	else:
		sys.exit('vcf file not found!')
	return vcf_df

def parse_split(chr_vcfs_dir_path, files_regex):
	printerr('parsing chromosome vcf files from ' + chr_vcfs_dir_path + '...')
	vcf_files = list(filter(re.compile(files_regex).match, os.listdir(chr_vcfs_dir_path)))

	vcf_list = []
	for vcf_file in vcf_files:
		
		printerr('parsing ' + vcf_file + '...')
		vcf_file_path = os.path.join(chr_vcfs_dir_path, vcf_file)
		index_length = pprogress.get_file_length(vcf_file_path)
		progress_index = pprogress.reset()

		with open(vcf_file_path, 'r') as inp:
			reader = vcf.VCFReader(inp)

			samples = reader.samples
			vcf_list.append(['variant_name'] + samples)

			for record in reader:
				variant_name = [vcfuid.rec_to_uid(record)]
				genotypes = []
				for s in range(len(samples)):
					gt = str_to_int(record.genotype(samples[s]).data[0])
					genotypes.append(gt)
				line = variant_name + genotypes

				vcf_list.append(line)
				progress_index = pprogress.pprogress(progress_index, index_length)
		
		vcf_array = np.array(vcf_list)
		vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
	return vcf_df