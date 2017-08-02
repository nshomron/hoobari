# --------- import modules ------------
# external
import re
import os
import sys
import subprocess
import requests
import sqlite3
import vcf, vcf.utils
import numpy as np
import pandas as pd
from json_commands import *
from pkl_commands import *
from collections import OrderedDict
# project's
import parse_gt
from stderr import *
import vcfuid
import pprogress
import position
import vcf_out
import preprocessing
<<<<<<< HEAD
from arguments import args
=======

# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')
parser.add_argument("-db", "--db_path", default = os.path.join(os.getcwd(), 'tmp_hb', 'hoobari.db'), help = 'path for vcf output')
parser.add_argument("-pkl", "--preprocessing_pkl_path", default = False, help = 'path to load preprocessing data from')
parser.add_argument("-model", "--model", default = 'simple', help = '	model for likelihoods calculation. possible values: "simple" \
									(Bayesian model based only on fetal fraction and parental genotypes), \
									"lengths" (use different fetal fraction per fragment length), \
									"origin" (use fragments that are very likely to be fetal, \
									based on other SNPs on these fragments)')

args = parser.parse_args()
>>>>>>> sqlite

# --------- pre-processing ----------
conn = sqlite3.connect(args.db_path)

# err_rate = 0.0003
# json_dir = os.path.join(args.tmp_dir, 'jsons')
# parents_gt = preprocessing.parse_parents_vcf(args.parents_vcf, fetal_sex = None, from_pkl = True, pkl_path = os.path.join(args.tmp_dir, 'parents_gt.pkl'))
# shared_fragments_dic, fetal_fragments_dic = preprocessing.create_fetal_and_shared_fragment_pools(	parents_gt,
# 													args.maternal_sample_name,
# 													args.paternal_sample_name,
# 													json_dir)
# total_fetal_fraction = preprocessing.calculate_total_fetal_fraction(shared_fragments_dic, fetal_fragments_dic)
# fetal_fractions_df = preprocessing.create_fetal_fraction_per_length_df(	shared_fragments_dic,
# 									fetal_fragments_dic,
# 									window = 3,
# 									max = 500,
# 									plot_dir = False)
# del shared_fragments_dic
# known_fetal_frags_dic = preprocessing.get_all_known_fetal_fragments(	parents_gt,
# 									args.maternal_sample_name,
# 									args.paternal_sample_name,
# 									fetal_fragments_dic,
# 									json_dir)

# pkl_paths = [os.path.join(args.tmp_dir, i) for i in ['err_rate', 'parents_gt', 'total_fetal_fraction', 'fetal_fractions_df', 'known_fetal_frags_dic']]
# pkl_save(err_rate, pkl_paths[0])
# pkl_save(parents_gt, pkl_paths[1])
# pkl_save(total_fetal_fraction, pkl_paths[2])
# pkl_save(fetal_fractions_df, pkl_paths[3])
# pkl_save(known_fetal_frags_dic, pkl_paths[4])

pkl_paths = [os.path.join(args.tmp_dir, i + '.pkl') for i in ['err_rate', 'parents_gt', 'total_fetal_fraction', 'fetal_fractions_df', 'known_fetal_frags_dic']]
err_rate = pkl_load(pkl_paths[0])
parents_gt = pkl_load(pkl_paths[1])
total_fetal_fraction = pkl_load(pkl_paths[2])
fetal_fractions_df = pkl_load(pkl_paths[3])
if args.model == 'origin':
	known_fetal_frags_dic = pkl_load(pkl_paths[4])
else:
	known_fetal_frags_dic = {}



cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)
if args.region:
	args_region_split = args.region.split(':')
	chrom = args_region_split[0]
	if len(args_region_split) > 1: 
		start, end = [int(i) for i in args_region_split[1].split('-')]
		cfdna_reader = cfdna_reader.fetch(chrom, start, end)
		parents_reader = parents_reader.fetch(chrom, start, end)
	else:
		cfdna_reader = cfdna_reader.fetch(chrom)
		parents_reader = parents_reader.fetch(chrom)

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
	if (tup[0] is None) or (tup[1] is None):
		print(tup)
	joint_probabilities = prediction = phred = probabilities_source = None

	cfdna_rec, parents_rec = tup
	
	if parents_rec:
		if cfdna_rec:
		
			# calculate priors
			maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
			paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
			priors = position.calculate_priors(maternal_gt, paternal_gt)
			

			if maternal_gt in (0,1,2):

				# calculate likelihoods
				likelihoods = position.calculate_likelihoods(	cfdna_rec,
										maternal_gt,
										args.tmp_dir,
										total_fetal_fraction,
										fetal_fractions_df,
										err_rate,
										known_fetal_frags_dic,
										args.model)

				# calculate posteriors
				joint_probabilities, prediction, phred, probabilities_source = position.calculate_posteriors(priors, likelihoods)
				
				if joint_probabilities is not None:
					# fetal information for the sample and FORMAT fields
					cfdna_geno_sample_dic = vcf_out.rec_sample_to_string(cfdna_rec, cfdna_id)
					if cfdna_geno_sample_dic != '.':
						cfdna_geno_sample_dic['GT'] = parse_gt.int_to_str(prediction)
						del cfdna_geno_sample_dic['GL']
						cfdna_geno_sample_dic['GJ'] = (','.join(str(round(p,2)) for p in list(joint_probabilities)))


					# parental information for INFO field
					parents_format = parents_rec.FORMAT
					
					if parents_rec.genotype(mother_id).data.GT != '.':
						matinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, mother_id).values()])
					if parents_rec.genotype(father_id).data.GT != '.':
						patinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, father_id).values()])
					else:
						matinfo = patinfo = '.'

					rec_info_dic = OrderedDict([	('PARENTS_FORMAT', parents_format),
									('MAT_INFO', matinfo),
									('PAT_INFO', patinfo),
									('PARENTS_QUAL', str(parents_rec.QUAL)),
									('PROB_SOURCE', probabilities_source)])



					# write var out (to file passed with -v or to output)
					vcf_out.print_var(cfdna_rec, phred, rec_info_dic, cfdna_geno_sample_dic, out_path = args.vcf_output)
				else:
					vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)
			else:
				vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)

	else:
		vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)
