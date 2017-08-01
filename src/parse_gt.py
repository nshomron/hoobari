# import modules
# external
import os
import re
import pandas as pd
import numpy as np
import vcf
# project's
from stderr import *
import vcfuid
import pprogress
from pkl_commands import *

def str_to_int(gt):
	if gt is not None and gt is not '.':
		gt_split = gt.split('/')
		gt_sum = int(gt_split[0]) + int(gt_split[1])
		if gt_sum in (0,1,2):
			return gt_sum
		else:
			return 'unsupported'
	else:
		return None

def int_to_str(gt):
	if gt == 0:
		string = '0/0'
	elif gt == 1:
		string = '0/1'
	elif gt == 2:
		string = '1/1'
	else:
		string = '.'
	return string

def parse(vcf_file_path, from_pkl = False, pkl_path = False):
	printerr('parsing ' + vcf_file_path + '...')
	if from_pkl and os.path.isfile(pkl_path):
		vcf_df = pkl_load(pkl_path)
	else:
		if os.path.isfile(vcf_file_path):
			index_length = pprogress.get_file_length(vcf_file_path)
			progress_index = pprogress.reset()
			percents = pprogress.percents(index_length)

			vcf_list = []
			reader = vcf.Reader(filename = vcf_file_path)

			samples = reader.samples
			vcf_list.append(['variant_name'] + samples + ['indel'])

			for record in reader:
				variant_name = [vcfuid.rec_to_uid(record)]
				genotypes = []
				for s in range(len(samples)):
					gt = str_to_int(record.genotype(samples[s]).data.GT)
					genotypes.append(gt)
				
				# flag: snp - 0 , indel - 1
				if record.INFO['TYPE'][0] in ('ins', 'del'):
					indel = 1
				else:
					indel = 0

				line = variant_name + genotypes + [indel]

				vcf_list.append(line)
				
				progress_index = pprogress.pprogress(progress_index, index_length, percents)
			
			vcf_array = np.array(vcf_list)
			vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
			pkl_save(vcf_df, pkl_path)
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

		reader = vcf.Reader(filename = vcf_file_path)

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