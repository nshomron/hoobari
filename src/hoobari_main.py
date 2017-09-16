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
from collections import OrderedDict
import pickle
import time
# project's
import parse_gt
from stderr import *
import vcfuid
import position
import vcf_out
import preprocessing
from arguments import args

sql_connection = sqlite3.connect(args.db)

# pre-processing
if os.path.isfile(args.preprocessing_pkl):
	time.sleep(5) # in case the file have just been created and is still being written
	with open(args.preprocessing_pkl, 'rb') as f:
		err_rate, total_fetal_fraction, fetal_fractions_df = pickle.load(f)
else:
	err_rate, total_fetal_fraction, fetal_fractions_df = preprocessing.run_full_preprocessing(	args.db,
													cores = args.cores,
													db_prefix = args.db_prefix,
													window = args.window,
													max_len = 500,
													plot = args.plot_lengths)
	if args.preprocessing_pkl:
		with open(args.preprocessing_pkl, 'wb') as f:
			pickle.dump((err_rate, total_fetal_fraction, fetal_fractions_df), f)

# create vcf files readers
cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)

# print header
input_command = ' '.join(sys.argv)
vcf_out.make_header(	cfdna_reader,
			parents_reader,
			input_command,
			args.fetal_sample_name,
			vcf_out.reserved_formats,
			output_path = args.vcf_output)

# fetch region
if args.region:
	region_split = re.split(':|-', args.region)
	chrom = region_split[0]
	start = int(region_split[1]) if len(region_split) > 1 else None
	end = int(region_split[2]) if len(region_split) > 2 else None
	try:
		cfdna_reader = cfdna_reader.fetch(chrom, start, end)
		parents_reader = parents_reader.fetch(chrom, start, end)
	except ValueError as e:
		errmessage = e.args[0]
		if errmessage.find('could not create iterator for region') == 0:
			sys.exit('warning! ' + errmessage + ', probably the input file does not contain any variants in the region.')


# get sample names from vcf files and from arguments
cfdna_id = cfdna_reader.samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name

# processing
co_reader = vcf.utils.walk_together(cfdna_reader, parents_reader)
for tup in co_reader:
#	if (tup[0] is None) or (tup[1] is None):
#		print(tup)
	joint_probabilities = prediction = phred = None

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
				posteriors, prediction, phred = position.calculate_posteriors(priors, likelihoods)
				normalized_likelihoods = position.likelihoods_to_phred_scale(likelihoods)

				# fetal information for the sample and FORMAT fields
				cfdna_geno_sample_dic = vcf_out.rec_sample_to_string(cfdna_rec, cfdna_id)
				if cfdna_geno_sample_dic != '.':
					cfdna_geno_sample_dic['GT'] = parse_gt.int_to_str(prediction)
					cfdna_geno_sample_dic['GL'] = (','.join(str(round(p,2)) for p in list(normalized_likelihoods)))
					cfdna_geno_sample_dic['PG'] = (','.join(str(round(p,5)) for p in list(priors)))
					cfdna_geno_sample_dic['PP'] = (','.join(str(round(p,5)) for p in list(posteriors)))




				# parental information for INFO field
				parents_format = parents_rec.FORMAT

				if parents_rec.genotype(mother_id).data.GT != '.':
					matinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, mother_id).values()])
				if parents_rec.genotype(father_id).data.GT != '.':
					patinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, father_id).values()])
				else:
					matinfo = patinfo = ':'.join(['.'] * (parents_rec.FORMAT.count(':') + 1))

				# INFO field of output vcf
				rec_info_list = []
				for s in (mother_id, father_id):
					for f in parents_rec.FORMAT.split(':'):
						if s == mother_id:
							pre = 'M'
						elif s == father_id:
							pre = 'P'
						parents_sample_field_data = parents_rec.genotype(s)[f]
						if type(parents_sample_field_data) is list:
							parents_sample_field_data = ','.join([str(i) for i in parents_sample_field_data])
						rec_info_list.append((pre + f, parents_sample_field_data))
				rec_info_list.append(('MPQ', str(parents_rec.QUAL)))
				rec_info_dic = OrderedDict(rec_info_list)


				# write var out (to file passed with -v or to output)
				vcf_out.print_var(cfdna_rec, phred, rec_info_dic, cfdna_geno_sample_dic, out_path = args.vcf_output)

			else:
				vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)

	else:
		vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)
