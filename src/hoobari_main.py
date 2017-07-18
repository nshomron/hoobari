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

class fetal_vcf():
	
	class header():

		def info_or_format(field_id, number, field_type, description, source=False):
			line_list = []
			line_list.append('##INFO=<ID=' + field_id)
			line_list.append('Number=' + number) # int, A, R, G, '.'
			line_list.append('Type=' + field_type)
			if source != False:
				line_list.append('Source="' + source + '"')
			line_list.append('Description="' + description + '">')

			return ','.join(line_list)
	
	def print_to_vcf(x, *args, **kargs):
		if args.vcf_output is not None:
			print(x, file = open(args.vcf_output, 'w'), *args, **kargs)
		else:
			print(x, *args, **kargs)


def make_vcf_header(cfdna_vcf_reader, parents_vcf_reader, fetal_sample_name, vcf_output_path):
	

	##fileformat=VCFv4.2
	##phasing=none

	if parents_vcf_reader.contigs == cfdna_vcf_reader.contigs:
		for f in parents_vcf_reader.contigs:


	for f in parents_vcf_reader.formats:
		write_fetal_vcf.header.info_or_format(	parents_vcf_reader.formats[f].id,
												parents_vcf_reader.formats[f].num,
												parents_vcf_reader.formats[f].type,
												parents_vcf_reader.formats[f].desc,
												'parental vcf file')



	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]



def write_var_to_vcf(variant_name, prediction, phred, pos_info_dic, format_and_gt_dic, vcf_output_path):

	vcf_list = []
	for row in posteriors_df.itertuples():

		row_list = []
		variant_name = row[0]
		
		# columns 1 - 5
		chrom, pos, ref, alt = parse_gt.uid_to_rec(variant_name)
		row_list += [chrom, pos, '.', ref, alt]
		
		# column 6-7
		row_list += [str(phred), '.']

		# column 8
		info_list = sorted([str(k) + '=' + str(info_dic[k]) for k in info_dic])
		row_list += [';'.join(info_list)]

		# columns 9-10
		format_list = []
		fetal_gt_list = []
		
		format_list.append('GT')
		fetal_gt_list.append(parse_gt.int_to_str(prediction))
		
		

		row_list += [':'.join(format_list)]
		row_list += [':'.join(fetal_gt_list)]

		vcf_list.append(row_list)
		
		print_var(variant_row)











#temp
err_rate = 0.0003

cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)
cfdna_id = cfdna_reader.samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name
co_reader = vcf.utils.walk_together(cfdna_reader, parents_reader)

for tup in co_reader:
	cfdna_rec, parents_rec = tup
	variant_name = vcfuid.rec_to_uid(cfdna_rec)
	
	# calculate priors
	maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	paternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	priors = position.calculate_priors(maternal_gt, paternal_gt)

	# calculate likelihoods
	likelihoods = position.calculate_likelihoods(variant, args.tmp_dir,	total_fetal_fraction, fetal_fractions_df, err_rate,	lengths = False, origin = False)

	# calculate posteriors
	posteriors, prediction, phred = position.calculate_posteriors(priors, likelihoods)

	info_dic = {'MATINFO_FORMAT': a,
				'MAT_INFO': B,
				'PAT_FORMAT': C,
				'PAT_INFO': D}

	gl = (','.join(str(p) for p in list(posteriors)))

	format_and_gt_dic = {	'GT': gt,
							'DP': dp,
							'AD': ad,
							'RO': ro,
							'QR': qr,
							'AO': ao,
							'QA': qa,
							'GL': gl}
