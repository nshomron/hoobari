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
from time import strftime
from json_commands import *
import argparse
from collections import OrderedDict
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

def print_info_or_format_row(field_id, number, field_type, description, source=False):
	line_list = []
	line_list.append('##INFO=<ID=' + field_id)
	line_list.append('Number=' + number) # int, A, R, G, '.'
	line_list.append('Type=' + field_type)
	if source:
		line_list.append('Source="' + source + '"')
	line_list.append('Description="' + description + '">')

	return ','.join(line_list)
	
def print_to_vcf(x, out_path = args.vcf_output, *args, **kargs):
	if out_path:
		print(x, file = open(out_path, 'w'), *args, **kargs)
	else:
		print(x, *args, **kargs)

def make_vcf_header(cfdna_vcf_reader, parents_vcf_reader, fetal_sample_name, infos, cfdna_geno_sample_dic):
	
	if cfdna_vcf_reader.metadata['reference'] != parents_vcf_reader.metadata['reference']:
		sys.exit('ERROR: the vcf files are not based on the same reference genome!')
	if cfdna_vcf_reader.contigs != parents_vcf_reader.contigs:
		printerr(	'Warning: cfdna and parental vcf files were aligned to same reference file \
					but their contigs are different')

	# print unique header fields
	print_to_vcf(	'##fileFormat=' + cfdna_vcf_reader.metadata['fileformat'],
					'##fileDate=' + strftime('%Y%m%d'),
					'##source=hoobari',
					'##phasing=none', # phasing is not yet supported
					'##reference=' + cfdna_vcf_reader.metadata['reference'],
					'##commandline="' + ' '.join(sys.argv) + '"',
					sep = '\n')
	
	# print contigs header fields
	cfdna_contigs_output = []
	for c in cfdna_vcf_reader.contigs.values():
		cfdna_contigs_output.append('##contig=<ID=' + str(c.id) + ',length=' + str(c.length) + '>')
	print_to_vcf('\n'.join(cfdna_contigs_output))

	# TODO ------------------
	# print info header fields
	for k in info_dic:
		print_info_or_format_row(	k,
									info_dic[k]['num'],
									info_dic[k]['type'],
									info_dic[k]['desc'],
									info_dic[k]['source'])

	# print format header fields
	cfdna_infos_names = [i[0] for i in cfdna_vcf_reader.formats.values()]
	for k in reserved_formats:
		if (k == 'GT') or (k == 'GJ'):
			k_source = '"hoobari"'
		else:
			k_source = '"cfdna vcf file"'

		print_info_or_format_row(	cfdna_vcf_reader.formats[k].id,
									vcf.Writer.counts[vcf_num_parents_vcf_reader.formats[k].num],
									cfdna_vcf_reader.formats[k].type,
									cfdna_vcf_reader.formats[k].desc,
									k_source)

	# ----------------------------

	# print column names
	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]
	print_to_vcf('\t'.join(vcf_columns))

def rec_sample_to_string(rec, sample):
	data = rec.genotype(sample).data
	
	gt = str(data.GT)
	dp = str(data.DP)
	ad = ','.join(str(i) for i in data.AD)
	ro = str(data.RO)
	qr = str(data.QR)
	ao = str(data.AO[0])
	qa = str(data.QA[0])
	gl = ','.join(str(i) for i in data.GL)

	format_and_gt_dic = OrderedDict([	('GT', gt),
										('DP', dp),
										('AD', ad),
										('RO', ro),
										('QR', qr),
										('AO', ao),
										('QA', qa),
										('GL', gl)])

	return format_and_gt_dic

def write_var_to_vcf(variant_name, prediction, phred, pos_info_dic, format_and_gt_dic):

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
		info_list = [str(k) + '=' + str(info_dic[k]) for k in info_dic]
		row_list += [';'.join(info_list)]

		# columns 9-10
		format_list = list(format_and_gt_dic.keys())
		fetal_gt_list = list(format_and_gt_dic.values())
		row_list += [':'.join(format_list)] + [':'.join(fetal_gt_list)]

		# merge all to one row string
		variant_row = '\t'.join(row_list)
		
		print_to_vcf(variant_row)

#temp
err_rate = 0.0003
infos = (, , , )
info_dic = OrderedDict([('MATINFO_FORMAT', {'num': '.', 'type': 'String', 'desc': '"Format of maternal sample ganotyping information"', 'source': '"parental vcf"'}),
						('MAT_INFO', {'num': '.', 'type': 'String', 'desc': '"Maternal sample ganotyping information"', 'source': '"parental vcf"'}),
						('PATINFO_FORMAT', {'num': '.', 'type': 'String', 'desc': '"Format of paternal sample ganotyping information"', 'source': '"parental vcf"'}),
						('PAT_INFO', {'num': '.', 'type': 'String', 'desc': '"Paternal sample ganotyping information"', 'source': '"parental vcf"'}),
						('PARENTS_QUAL', {'num': '.', 'type': 'Float', 'desc': '"Parental samples QUAL score"', 'source': '"parental vcf"'})])
reserved_formats = ('GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL')

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
	maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
	paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
	priors = position.calculate_priors(maternal_gt, paternal_gt)

	# calculate likelihoods
	likelihoods = position.calculate_likelihoods(variant, args.tmp_dir,	total_fetal_fraction, fetal_fractions_df, err_rate,	lengths = False, origin = False)

	# calculate posteriors
	joint_probabilities, prediction, phred = position.calculate_posteriors(priors, likelihoods)

	# parental information for INFO field
	parents_format = parents_reader.FORMAT
	matinfo = ':'.join(str(i) for i in rec_sample_to_string(parents_rec, mother_id).values())
	patinfo = ':'.join(str(i) for i in rec_sample_to_string(parents_rec, father_id).values())
	rec_info_dic = OrderedDict([('MATINFO_FORMAT', parents_format),
								('MAT_INFO', matinfo),
								('PATINFO_FORMAT', parents_format),
								('PAT_INFO', patinfo),
								('PARENTS_QUAL', str(parents_rec.QUAL))])

	# fetal information for the sample and FORMAT fields
	cfdna_geno_sample_dic = rec_sample_to_string(cfdna_rec, cfdna_id)
	cfdna_geno_sample_dic['GT'] = parse_gt.int_to_str(prediction)
	cfdna_geno_sample_dic['GJ'] = (','.join(str(p) for p in list(joint_probabilities)))
	del cfdna_geno_sample_dic['GL']

	# write var out (to file passed with -v or to output)
	write_var_to_vcf(variant_name, prediction, phred, rec_info_dic, cfdna_geno_sample_dic)