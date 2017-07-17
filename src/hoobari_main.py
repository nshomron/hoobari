# --------- import modules ------------
# external
import re
import os
import sys
import subprocess
import requests
import vcf, vcf.utils
import numpy as np
import pandas as pd
from json_commands import *
import argparse
# project's
import parse_gt
from stderr import printerr
import vcfuid
import pprogress
import position

# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-parents", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-tmp", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-v", "--vcf_output", help = 'path for vcf output')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')

args = parser.parse_args()

printerr(args)

# --------- functions ----------








def parse_vcf_to_array(vcf_file_path, pkl_file_path):
	stderr('parsing ' + vcf_file_path + '...')
	if os.path.isfile(pkl_file_path):
		vcf_df = pkl_load(pkl_file_path)
	else:
		if os.path.isfile(vcf_file_path):
			index_length = int(subprocess.getoutput(' '.join(['grep','-v','^#',vcf_file_path,'|','wc','-l'])).split()[0])
			progress_index = 0

			vcf_list = []
			with open(vcf_file_path, 'r') as inp:
				reader = vcf.VCFReader(inp)

				samples = reader.samples
				vcf_list.append(['variant_name'] + samples)

				for record in reader:
					variant_name = [vcf_rec_to_var_uid(record)]
					genotypes = []
					for s in range(len(samples)):
						gt = gt_string_to_int(record.genotype(samples[s]).data[0])
						genotypes.append(gt)
					line = variant_name + genotypes

					vcf_list.append(line)
					progress_index = print_progress(progress_index, index_length)
				
				vcf_array = np.array(vcf_list)
				vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
			pkl_save(vcf_df, pkl_file_path)
		else:
			sys.exit('vcf file not found!')
	return vcf_df


print_var(variant_row, *args, **kargs):
	if args.vcf_output is not None:
		print(variant_row, file = open(args.vcf_output, 'w'), *args, **kargs)
	else:
		print(variant_row, *args, **kargs)

def make_vcf_header(fetal_sample_name, vcf_output_path):
	
	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]

def write_var_to_vcf(vcf_output_path):

	posteriors_df['gt'] = posteriors_df.ix[:,0:3].idxmax(axis=1)
	posteriors_df['phred'] = np.where(posteriors_df['gt'] == 0, (-10) * np.log10(1 - posteriors_df[0]), (-10) * np.log10(posteriors_df[0]))

	vcf_list = []
	for row in posteriors_df.itertuples():

		row_list = []
		variant_name = row[0]
		
		# columns 1 - 5
		var_split = re.split(r':|_|/', variant_name)
		row_list += var_split[0:2] + ['.'] + var_split[2:4]
		
		# column 6-7
		row_list += [str('phred')] + ['.']

		# column 8
		info_list = []
		info_list.append('MATINFO_FORMAT=' + '.')
		info_list.append('MAT_INFO=' + '.')
		info_list.append('PAT_FORMAT=' + '.')
		info_list.append('PAT_INFO=' + '.')
		info_list.append('POSTERIORS=' + (','.join(str(p) for p in list(posteriors))))
		row_list += [';'.join(info_list)]

		# columns 9-10
		format_list = []
		fetal_gt_list = []
		
			# GT:DP:AD:RO:QR:AO:QA:GL

		format_list.append('GT')
		fetal_gt_list.append(gt_int_to_string(posteriors_df.ix[variant_name, 'gt']))
		
		row_list += [':'.join(format_list)]
		row_list += [':'.join(fetal_gt_list)]

		vcf_list.append(row_list)
		
		print_var(variant_row)











#temp
err_rate = 0.0003

readers = [vcf.Reader(open(args.cfdna_vcf, 'rb')), vcf.Reader(open(args.parents_vcf, 'rb'))]
cfdna_id = readers[0].samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name
reader = vcf.utils.walk_together(readers)

for tup in reader:
	cfdna_rec, parents_rec = tup
	variant_name = vcfuid.rec_to_uid(cfdna_rec)
	
	# calculate priors
	maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	paternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	priors = position.calculate_priors(maternal_gt, paternal_gt)

	# calculate likelihoods
	likelihoods = position.calculate_likelihoods(variant, args.tmp_dir,	total_fetal_fraction, fetal_fractions_df, err_rate,	lengths = False, origin = False)

	# calculate posteriors
	posteriors, phred = position.calculate_posteriors(priors, likelihoods)

