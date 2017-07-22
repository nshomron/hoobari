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
import vcf_out

# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')


args = parser.parse_args()

printerr(args)

# --------- functions ----------


err_rate = 0.0003
cfdna_reader = vcf.Reader(filename = args.cfdna_vcf)
parents_reader = vcf.Reader(filename = args.parents_vcf)
cfdna_id = cfdna_reader.samples[0]
mother_id = args.maternal_sample_name
father_id = args.paternal_sample_name

vcf_out.make_header(	cfdna_reader,
			parents_reader,
			args.fetal_sample_name,
			vcf_out.info_dic,
			vcf_out.reserved_formats,
			out_path = args.vcf_output)

co_reader = vcf.utils.walk_together(cfdna_reader, parents_reader)
for tup in co_reader:
	cfdna_rec, parents_rec = tup
	variant_name = vcfuid.rec_to_uid(cfdna_rec)
	
	# calculate priors
	maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data.GT)
	paternal_gt = parse_gt.str_to_int(parents_rec.genotype(father_id).data.GT)
	priors = position.calculate_priors(maternal_gt, paternal_gt)

	# calculate likelihoods
	likelihoods = position.calculate_likelihoods(variant, args.tmp_dir, total_fetal_fraction, fetal_fractions_df, err_rate,	lengths = False, origin = False)

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
	cfdna_geno_sample_dic['GJ'] = (','.join(str(p) for p in list(joint_probabilities)))
	del cfdna_geno_sample_dic['GL']

	# write var out (to file passed with -v or to output)
	vcf_out.print_var(variant_name, prediction, phred, rec_info_dic, cfdna_geno_sample_dic, out_path = args.vcf_output)