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
import preprocessing_lowmem as preprocessing
# import preprocessing
from arguments import args

sql_connection = sqlite3.connect(args.db)

# pre-processing
err_rate, total_fetal_fraction, fetal_fractions_df = preprocessing.run_full_preprocessing(args.db, cores = args.cores, db_prefix = args.db_prefix, window = 3, max_len = 500, plot = args.plot_lengths)

# processing
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
#	if (tup[0] is None) or (tup[1] is None):
#		print(tup)
	joint_probabilities = prediction = phred = probabilities_source = None

	cfdna_rec, parents_rec = tup
	
	printverbose('parents_rec: ', parents_rec)
	if parents_rec:
		printverbose('cfdna_rec: ', cfdna_rec)
		if cfdna_rec:
		
			# calculate priors
			maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
			paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
			priors = position.calculate_priors(maternal_gt, paternal_gt)
			
			printverbose(maternal_gt, paternal_gt)
			if maternal_gt in (0,1,2):

				# calculate likelihoods
				likelihoods = position.calculate_likelihoods(	cfdna_rec,
										maternal_gt,
										total_fetal_fraction,
										fetal_fractions_df,
										err_rate,
										sql_connection,
										args.model)

				# calculate posteriors
				joint_probabilities, prediction, phred, probabilities_source = position.calculate_posteriors(priors, likelihoods)
				
				printverbose('joint_probabilities', joint_probabilities)
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
