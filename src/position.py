# --------- import modules ------------
# external
import os
import requests
import numpy as np
import time
import argparse
from multiprocessing import Pool, cpu_count
import ctypes
from decimal import *
# project's
import parse_gt
from stderr import *
import vcfuid

# --------- global -------------
de_novo = 1.2e-8
valid_gts = (0,1,2)
default_decimal_prec = getcontext().prec

# --------- functions ----------
def calculate_priors(maternal_gt, paternal_gt):

	'''
	Calculate prior probabilities for each of the 3 possible fetal genotypes (1/1, 0/1 and 0/0), given the parents genotypes.
	Mother: Calculations below assume maternal genotype = '0/1', since these are the positions of interest. Therefore the
	probability to inherit the alternate allele from her is 0.5.
	Father: If the father is 1/1 (2) the probability to inherit the alternate allele is 1. if he's 0/1 (1) it's 0.5,
	and if he's 0/0 (0) it's 0. therefore it's always the (sum_of_alternate_alleles/2), which is marked is f.
	P(fetus = 1/1) = p_inherit_alt_from_mother * p_inherit_alt_from_father = 0.5 * f
	P(fetus = 0/1) = (p_inherit_alt_from_mother * p_inherit_ref_from_father) + (p_inherit_ref_from_mother * p_inherit_alt_from_father) =
	0.5 * (1 - f) + 0.5 * f = 0.5 * (1 - f + f) = 0.5
	P(fetus = 0/0) = p_inherit_ref_from_mother * p_inherit_ref_from_father = 0.5 * (1 - f)
	==>
	[0.5f, 0.5, 0.5(1-f)]
	'''

	# TODO: de-novo mutation rate per site
	# TODO: de-novo mutation rate per
	# TODO: add error rate to prior

	if (maternal_gt in valid_gts) and (paternal_gt in valid_gts):
		p_maternal_alt = maternal_gt / 2
		p_paternal_alt = paternal_gt / 2
		priors = np.array([	(1-p_maternal_alt)*(1-p_paternal_alt),
					p_maternal_alt*(1-p_paternal_alt) + (1-p_maternal_alt)*p_paternal_alt,
					p_maternal_alt*p_paternal_alt])
	else:
		priors = np.array([1/3, 1/3, 1/3])

	printverbose('priors:', priors)

	return priors

def sqname(bam_rec):
	return(bam_rec.query_name + ':' + str(int(not bam_rec.is_reverse)))


def fetch_bam_records(bam_reader, chrom, position):
	start_pos = max(0, int(position) - 1000)
	chrom = 'chr' + str(chrom)
	end_pos = min(bam_reader.lengths[bam_reader.references.index(chrom)], int(position) + 1000)
	bam_records_at_position = bam_reader.fetch(chrom, start_pos, end_pos) # include a flanking region, since there's local realignment
	rec_at_position_dic = {}
	for rec in bam_records_at_position:
		# rec_at_position_dic[sqname(rec)] = rec
		if rec.query_name in rec_at_position_dic:
			rec_at_position_dic[rec.query_name].append(rec)
		else:
			rec_at_position_dic[rec.query_name] = [rec]
	# print(rec_at_position_dic)
	return rec_at_position_dic

def fetch_mate(bam_reader, bam_rec):
	recs_at_mate_pos = bam_reader.fetch(	bam_rec.next_reference_name,
						bam_rec.next_reference_start,
						bam_rec.next_reference_start + 1)
	for rec in recs_at_mate_pos:
		if rec.query_name == bam_rec.query_name and rec.reference_start != bam_rec.reference_start:
			return(rec)

	return(None)
	
def calculate_fragment_i(frag_genotype, maternal_gt, ref, alt, p_fet, err_rate):

	if frag_genotype == alt:
		p_maternal_alt = max(err_rate, maternal_gt / 2)
		frag_i_likelihoods = np.array([	0*p_fet + p_maternal_alt*(1-p_fet), # fetal is 0/0
						0.5*p_fet + p_maternal_alt*(1-p_fet), # fetal is 0/1
						1*p_fet + p_maternal_alt*(1-p_fet)]) # fetal is 1/1

	elif frag_genotype == ref:
		p_maternal_ref = max(err_rate, 1 - (maternal_gt / 2))
		frag_i_likelihoods = np.array([	1*p_fet + p_maternal_ref*(1-p_fet), # fetal is 0/0
						0.5*p_fet + p_maternal_ref*(1-p_fet), # fetal is 0/1
						0*p_fet + p_maternal_ref*(1-p_fet)]) # fetal is 1/1

	else:
		frag_i_likelihoods = np.array([1, 1, 1])

	return frag_i_likelihoods

def calculate_likelihoods(
	rec,
	cfdna_bam_reader,
	maternal_gt,
	total_fetal_fraction,
	fetal_fractions_df,
	err_rate,
	vardb,
	model,
	max_depth,
	**kwargs):

	'''
	calculate likelihoods for each of the 3 possible fetal genotype, based on the prbability that
	fragment_i at the position, will show a certain allele (ref or alt), given other factors of
	the model, such as the maternal genotype, fragment length and the fetal genotype (which is unknown,
	so we check for all possibilities - 1/1, 0/1 and 0/0)
	'''

	chrom, pos, ref, alt = rec.CHROM.replace('chr', ''), str(rec.POS), rec.REF, str(rec.ALT[0])

	variant_len = len(ref) - len(alt)

	printverbose(chrom, pos)
	pos_data = vardb.getPositionVariants(chrom, pos, model).fetchall()
	printverbose('position_data:')
	printverbose(pos_data)

	first = 1
	if (len(pos_data) > 0) and (maternal_gt in valid_gts):

		for frag_genotype, frag_length, frag_is_fetal in pos_data:

			# if variant_len > 0: # deletion ####
			frag_length = max(int(frag_length) - variant_len, 0)

			# get fetal fraction
			if (model == 'origin') and (frag_is_fetal == 1):
				printverbose('fragment was marked as fetal by another variant')
				p_fetal = 1 - err_rate
			elif (model in ('lengths', 'origin')) and (frag_length in fetal_fractions_df.index.values):
				p_fetal = fetal_fractions_df[frag_length]
			else:
				p_fetal = total_fetal_fraction

			frag_i_likelihoods = calculate_fragment_i(frag_genotype,
								maternal_gt,
								ref,
								alt,
								p_fetal,
								err_rate)

			#frag_i_likelihoods[frag_i_likelihoods == 0] = 0.003
			frag_i_likelihoods = np.log(frag_i_likelihoods, dtype = np.float64)
			if first:
				sum_log_fragments_likelihoods_df = frag_i_likelihoods
				first = 0
			else:
				sum_log_fragments_likelihoods_df = np.add(sum_log_fragments_likelihoods_df, frag_i_likelihoods)
			printverbose(p_fetal, frag_i_likelihoods, sum_log_fragments_likelihoods_df, sep = '\t')

	else:
		sum_log_fragments_likelihoods_df = np.log([1, 1, 1], dtype = np.float64)
	printverbose(sum_log_fragments_likelihoods_df)

	return sum_log_fragments_likelihoods_df

def calculate_phred(joint_probabilities):
	libphred = np.ctypeslib.load_library('libphred.so','/groups/nshomron/tomr/projects/cffdna/hoobari/src/') #TODO - fix to relative path
	libphred.calculatePhred.restype = ctypes.c_longdouble
	cdubs = joint_probabilities.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble))

	return libphred.calculatePhred(cdubs)

def simple_qual_calculation(posteriors):

	if posteriors[0] == 0:
		return '1e+06'
	else:
		qual = -10 * np.log10(posteriors[0])

		if qual == 0:
			return int(0)
		elif qual < 1:
			return '{:0.3e}'.format(qual)
		else:
			return round(qual,2)

def likelihoods_to_phred_scale(likelihoods):

	log10_likelihoods = likelihoods / np.log(10)
	normalized_log10_likelihoods = log10_likelihoods - np.max(log10_likelihoods)
	phred_scaled_normalized_likelihoods = -10 * (normalized_log10_likelihoods)
	phred_scaled_normalized_likelihoods[phred_scaled_normalized_likelihoods == -0.0] = 0
	printverbose('normalized_likelihoods', phred_scaled_normalized_likelihoods)

	return phred_scaled_normalized_likelihoods

def calculate_posteriors(var_priors, var_likelihoods):

	printverbose(var_priors, var_likelihoods, sep = '\n')
	# sum priors and likelihoods
	# if there are no priors (for instance if parental genotypes at positions are missing),
	# take only likelihoods


	# parents genotypes might give a prediction but if it's not supported by cfdna it's only indirect - good or not? not for de-novo...
	# it is good for recessive disease!

	var_priors = np.log(var_priors)
	printverbose('log(priors):', var_priors)

	joint_probabilities = np.add(var_priors, var_likelihoods, dtype = np.float64)
	prediction = joint_probabilities.argmax()

	joint_probabilities_normalized = joint_probabilities - np.min(joint_probabilities[np.isfinite(joint_probabilities)])
	exp_joint_probabilities = np.exp(joint_probabilities_normalized)
	printverbose('joint_probabilities_normalized:', joint_probabilities_normalized)
	posteriors = np.asarray((exp_joint_probabilities / np.sum(exp_joint_probabilities)), dtype = np.float64)
	printverbose('posteriors:', posteriors)
	if np.inf in exp_joint_probabilities:
		getcontext().prec = default_decimal_prec
		start_time = time.time()
		while 1 or np.nan in posteriors:
			if getcontext().prec <= 2**11:
				printverbose('precision:', str(getcontext().prec))
				exp_joint_probabilities = np.array([Decimal(i).exp() for i in joint_probabilities_normalized])
				posteriors = np.array(exp_joint_probabilities / np.sum(exp_joint_probabilities))
				printverbose('posteriors:', posteriors)
				getcontext().prec *= 2
			else:
				break
		printverbose(time.time() - start_time)

	qual = simple_qual_calculation(posteriors)

	printverbose('Posteriors:', posteriors)
	printverbose('Predicted genotype:', prediction)
	printverbose('QUAL:', qual)

	return (posteriors, prediction, qual)
