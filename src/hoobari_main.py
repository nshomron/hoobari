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

# connect to the database that was created during the first analysis of the cfDNA sample
sql_connection = sqlite3.connect(args.db)

# pre-processing
# calculate the total fetal fraction and a table of fetal-fraction per fragment size
# calcuation results can be saved to a pkl file, and read from a pkl as well
# TODO: note that the error rate isn't actually used in the model yet
if os.path.isfile(args.preprocessing_pkl):
	time.sleep(5) # in case the file pkl have just been created and is still being written
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

# create vcf files iterators
cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)

# print header of output vcf file
input_command = ' '.join(sys.argv)
vcf_out.make_header(	cfdna_reader,
			parents_reader,
			input_command,
			args.fetal_sample_name,
			vcf_out.reserved_formats,
			output_path = args.vcf_output)

# fetch region, if a region was specified
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
		if 'could not create iterator for region' in errmessage:
			sys.exit('warning! ' + errmessage + ', probably the input file does not contain any variants in the region.')


# get sample of cfdna from its vcf file, and of the parents from the input arguments
cfdna_id = cfdna_reader.samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name

# processing positions
# iterate on both vcf files and return a tuple for each position, that contains its record from each vcf file.
# if there is no vcf record for a certain position in one of the files, a None will appear in the tuple instead.
co_reader = vcf.utils.walk_together(cfdna_reader, parents_reader)
for tup in co_reader:
	# reset prediction and QUAL
	prediction = qual = None

	cfdna_rec, parents_rec = tup
	
	printverbose('parents_rec: ', parents_rec)
	printverbose('cfdna_rec: ', cfdna_rec)
	
	if parents_rec and cfdna_rec:
		
		# calculate priors for the position
		maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
		paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
		priors = position.calculate_priors(maternal_gt, paternal_gt)
		
		printverbose(maternal_gt, paternal_gt)
		
		# for now, only positions where the mother is 0/0, 0/1 or 1/1 are supported
		if maternal_gt in (0,1,2):

			# calculate likelihoods for the position
			likelihoods = position.calculate_likelihoods(	cfdna_rec,
									maternal_gt,
									total_fetal_fraction,
									fetal_fractions_df,
									err_rate,
									sql_connection,
									args.model)

			# calculate posteriors for the position
			posteriors, prediction, qual = position.calculate_posteriors(priors, likelihoods)
			
			
			## process the output entry
			
			# create normalized likelihoods for the output vcf
			normalized_likelihoods = position.likelihoods_to_phred_scale(likelihoods)

			# fetal information for the sample and FORMAT fields
			cfdna_geno_sample_dic = vcf_out.rec_sample_to_string(cfdna_rec, cfdna_id)
			if cfdna_geno_sample_dic != '.':
				cfdna_geno_sample_dic['GT'] = parse_gt.int_to_str(prediction)
				cfdna_geno_sample_dic['GL'] = (','.join(str(round(p,2)) for p in list(normalized_likelihoods)))
				cfdna_geno_sample_dic['PG'] = (','.join(str(round(p,5)) for p in list(priors)))
				cfdna_geno_sample_dic['PP'] = (','.join(str(round(p,5)) for p in list(posteriors)))	


			# TODO: MARKED FOR DELETION
			# # process the parental information for the INFO field
			# parents_format = parents_rec.FORMAT
			# parents_format_len = parents_rec.FORMAT.count(':') + 1
			# empty_info = ':'.join(['.'] * parents_format_len)

			# # get maternal info
			# if parents_rec.genotype(mother_id).data.GT != '.':
			# 	matinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, mother_id).values()])
			# else:
			# 	matinfo = empty_info
			
			# # get paternal info
			# if parents_rec.genotype(father_id).data.GT != '.':
			# 	patinfo = ':'.join([str(i) for i in vcf_out.rec_sample_to_string(parents_rec, father_id).values()])
			# else:
			# 	patinfo = empty_info


			# for each parent, for all the data in its sample, create an instance that will be printed in the output INFO
			rec_info_list = []
			for parent_sample in (mother_id, father_id):
				for field in parents_rec.FORMAT.split(':'):
					if parent_sample == mother_id:
						prefix = 'M'
					elif parent_sample == father_id:
						prefix = 'P'
					parents_sample_field_data = parents_rec.genotype(parent_sample)[field]
					if type(parents_sample_field_data) is list: # some fields contain a few values
						parents_sample_field_data = ','.join([str(i) for i in parents_sample_field_data])
					rec_info_list.append((prefix + field, parents_sample_field_data))
			rec_info_list.append(('MPQ', str(parents_rec.QUAL)))
			rec_info_dic = OrderedDict(rec_info_list)

			# write var out (to file passed with -v or to output)
			vcf_out.print_var(cfdna_rec, qual, rec_info_dic, cfdna_geno_sample_dic, out_path = args.vcf_output)

		else:
			vcf_out.unsupported_position(cfdna_rec, out_path = args.vcf_output)

	else:
		if parents_rec is None:
			empty_rec = cfdna_rec
		else:
			empty_rec = parents_rec
		vcf_out.unsupported_position(empty_rec, out_path = args.vcf_output)