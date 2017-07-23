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
from pkl_commands import *
import argparse
from collections import OrderedDict
# project's
import parse_gt
from stderr import printerr
import vcfuid
import pprogress
import position
import vcf_out
import preprocessing

# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')
parser.add_argument("-pkl", "--preprocessing_pkl_path", default = False, help = 'path to load preprocessing data from')
parser.add_argument("-model", "--model", default = 'simple', help = '	model for likelihoods calculation. possible values: "simple" \
									(Bayesian model based only on fetal fraction and parental genotypes), \
									"lengths" (use different fetal fraction per fragment length), \
									"origin" (use fragments that are very likely to be fetal, \
									based on other SNPs on these fragments)')

args = parser.parse_args()

# --------- pre-processing ----------

if not args.preprocessing_pkl_path:
	err_rate = 0.0003
	parents_gt = preprocessing.parse_parents_vcf(args.parents_vcf, fetal_sex = None)
	shared_fragments_dic, fetal_fragments_dic = preprocessing.create_fetal_and_shared_fragment_pools(	parents_gt,
														args.maternal_sample_name,
														args.paternal_sample_name)
	total_fetal_fraction = preprocessing.calculate_total_fetal_fraction(shared_fragments_dic, fetal_fragments_dic)
	known_fetal_frags_dic = preprocessing.get_all_known_fetal_fragments(parents_gt, fetal_fragments_dic)
	fetal_fractions_df = preprocessing.create_fetal_fraction_per_length_df(	shared_fragments_dic,
										fetal_fragments_dic,
										window = 3,
										max = 500,
										plot_dir = False)

	pkl_save([err_rate, parents_gt, total_fetal_fraction, known_fetal_frags_dic, fetal_fractions_df], os.path.join(args.tmp_dir, 'pre_processing.pkl'))
else:
	err_rate, parents_gt, total_fetal_fraction, known_fetal_frags_dic, fetal_fractions_df = pkl_load(args.use_preprocessing_from_pkl)

cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)
cfdna_id = cfdna_reader.samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name

input_command = ' '.join(sys.argv) + '"'
vcf_out.make_header(	cfdna_reader,
			parents_reader,
			input_command,
			args.fetal_sample_name,
			vcf_out.info_dic,
			vcf_out.reserved_formats,
			output_path = args.vcf_output)

co_reader = vcf.utils.walk_together(cfdna_reader, parents_reader)
for tup in co_reader:
	cfdna_rec, parents_rec = tup
	
	if parents_rec:
		if cfdna_rec:
			variant_name = vcfuid.rec_to_uid(cfdna_rec)
			
			# calculate priors
			maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
			paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
			priors = position.calculate_priors(maternal_gt, paternal_gt)

			# calculate likelihoods
			likelihoods = position.calculate_likelihoods(	cfdna_rec,
									args.tmp_dir,
									total_fetal_fraction,
									fetal_fractions_df,
									err_rate,
									known_fetal_frags_dic,
									args.model)

			# calculate posteriors
			joint_probabilities, prediction, phred = position.calculate_posteriors(priors, likelihoods)

			# parental information for INFO field
			parents_format = parents_reader.FORMAT
			matinfo = ':'.join(str(i) for i in rec_sample_to_string(parents_rec, mother_id).values())
			patinfo = ':'.join(str(i) for i in rec_sample_to_string(parents_rec, father_id).values())
			rec_info_dic = OrderedDict([	('MATINFO_FORMAT', parents_format),
							('MAT_INFO', matinfo),
							('PATINFO_FORMAT', parents_format),
							('PAT_INFO', patinfo),
							('PARENTS_QUAL', str(parents_rec.QUAL))])

			# fetal information for the sample and FORMAT fields
			cfdna_geno_sample_dic = rec_sample_to_string(cfdna_rec, cfdna_id)
			cfdna_geno_sample_dic['GT'] = parse_gt.int_to_str(prediction)
			del cfdna_geno_sample_dic['GL']
			cfdna_geno_sample_dic['GJ'] = (','.join(str(p) for p in list(joint_probabilities)))

			# write var out (to file passed with -v or to output)
			vcf_out.print_var(cfdna_rec, phred, rec_info_dic, cfdna_geno_sample_dic, out_path = args.vcf_output)
	
	else:
		vcf_row = [	cfdna_rec.CHROM,
				str(cfdna_rec.POS),
				'.',
				cfdna_rec.REF,
				cfdna_rec.ALT,
				'.',
				'.',
				'MATINFO_FORMAT=.;MAT_INFO=.;PATINFO_FORMAT=.;PAT_INFO=.;PARENTS_QUAL=.',
				'.',
				'.']
		vcf_out.printvcf(variant_row, out_path = args.vcf_output)