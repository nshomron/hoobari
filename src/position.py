# --------- import modules ------------
# external
import re
import os
import sys
import subprocess
import requests
import vcf
import numpy as np
import pandas as pd
from json_commands import *
import argparse
from multiprocessing import Pool, cpu_count
import ctypes
# project's
import parse_gt
from stderr import *
import vcfuid

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
	if (maternal_gt in (0,1,2)) and (paternal_gt in (0,1,2)):
		p_maternal_alt = maternal_gt / 2
		p_paternal_alt = paternal_gt / 2
		priors_source = 'parents_vcf'
		priors = [	(1-p_maternal_alt)*(1-p_paternal_alt),
				p_maternal_alt*(1-p_paternal_alt) + (1-p_maternal_alt)*p_paternal_alt,
				p_maternal_alt*p_paternal_alt]

		for i in range(len(priors)):
			if priors[i] == 0:
				priors[i] = np.log(0.000000012) # TODO: None or some very small value/s?
			else:
				priors[i] = np.log(priors[i])

	# elif not maternal_gt and not paternal_gt:
	# 	priors = [0.25, 0.5, 0.25]
	# 	priors_source = 'naive'		
	# elif not maternal_gt:
	# 	p_paternal_alt = paternal_gt / 2
	# 	p_maternal_alt = 0.5
	# 	priors_source = 'only_paternal'
	# elif not paternal_gt:
	# 	p_maternal_alt = maternal_gt / 2
	# 	priors_source = 'only_maternal'
	elif (maternal_gt == 'unsupported') or (paternal_gt == 'unsupported'):
		priors = 'more than one alternate allele, not yet supported'
		priors_source = 'unsupported'
	else:
		priors = [None, None, None]
		priors_source = 'no_priors'
	# if priors_source != 'naive':

	# else:
	# 	get maf or af or ldaf
	# 	maf = retrieve_maf_from_ensembl(variant_name)
	# 	if maf[1] is not None:
	# 		p_maternal_alt = p_paternal_alt = maf[1]
	# 		priors_source = maf[0]




	return (priors, priors_source)

def calculate_fragment_i(frag_genotype, maternal_gt, ref, alt, f, err_rate):
	'''
	probabilities of each fragment i to show the reference allele, given one of the 3 possible fetal genotypes.
	P(frag_i | fetal_genotype) = P(frag_i | frag from fetus)P(frag from fetus | fetal_genotype) + P(frag_i | frag from mother)P(frag from mother | fetal_genotype)
	the first expression is the fetal quantity of the allele, and the second is the maternal quantity of the allele.
	example: fetal fraction (ff is 0.1, fetus is aa, mother Aa. Therefore the fetus donates 0.1 fragments with genotype a, and half of the maternal fragments,
	which is (1-0.1)/2 = 0.45, also donate genotype a => ff + (1-ff)*0.5 = 2ff*0.5 + (1-ff)*0.5 = 0.5(1 - ff + 2ff) = 0.5(1 + ff)
	if fetus is aA - 0.5*ff + 0.5*(1-ff) = 0.5(ff+1-ff) = 0.5
	if fetus is AA - 0*ff + 0.5(1-ff) = 0.5(1-ff)
	'''

	# fetal genotypes: 0,1,2
	# print(frag_genotype)
	if frag_genotype == alt:
		p_maternal_alt = maternal_gt / 2
		frag_i_likelihoods = [	0*f + p_maternal_alt*(1-f), # fetal is 0/0
					0.5*f + p_maternal_alt*(1-f), # fetal is 0/1
					1*f + p_maternal_alt*(1-f)] # fetal is 1/1
	elif frag_genotype == ref:
		p_maternal_ref = 1 - (maternal_gt / 2)
		frag_i_likelihoods = [	1*f + p_maternal_ref*(1-f), # fetal is 0/0
					0.5*f + p_maternal_ref*(1-f), # fetal is 0/1
					0*f + p_maternal_ref*(1-f)] # fetal is 1/1
	else:
		frag_i_likelihoods = None
	
	return frag_i_likelihoods

def calculate_likelihoods(
	rec,
	maternal_gt,
	tmp_dir,
	total_fetal_fraction,
	fetal_fractions_df,
	err_rate,
	known_fetal_qnames_dic,
	model,
	**kwargs):

	'''
	calculate likelihoods for each of the 3 possible fetal genotype, based on the prbability that
	fragment_i at the position, will show a certain allele (ref or alt), given other factors of
	the model, such as the maternal genotype, fragment length and the fetal genotype (which is unknown,
	so we check for all possibilities - 1/1, 0/1 and 0/0)
	'''
	sum_log_fragments_likelihoods_df = [None, None, None]

	chrom, pos, ref, alt = rec.CHROM, str(rec.POS), rec.REF, str(rec.ALT[0])
	
	variant_len = len(ref) - len(alt)

	snp_json_path = os.path.join(tmp_dir, 'jsons', chrom, pos + '.json')
	pos_data = pd.DataFrame(json_load(snp_json_path))

	if (pos_data is not None) and (maternal_gt in (0,1,2)):

		fragments_likelihoods_list = []
		for row in pos_data.itertuples():
			frag_genotype = row[1]
			frag_length = max(int(row[2]) - variant_len, 0)
			frag_qname = row[3]
			
			# get fetal fraction
			if (model == 'origin') and (frag_qname in known_fetal_qnames_dic):
				ff = 0.7
			elif (model in ('lengths', 'origin')) and (frag_length in fetal_fractions_df.index.values):
				ff = fetal_fractions_df[frag_length]
			else:
				ff = total_fetal_fraction

			frag_i_likelihood_list = calculate_fragment_i(frag_genotype, maternal_gt, ref, alt, ff, err_rate)
			if frag_i_likelihood_list is not None:
				for i in range(len(frag_i_likelihood_list)):
					if frag_i_likelihood_list[i] == 0: # 0 would cause -inf after log, and the sum would also be -inf
						frag_i_likelihood_list[i] = 0.03 # TODO: change this when error is inserted to the model
				fragments_likelihoods_list.append(frag_i_likelihood_list)

		if len(fragments_likelihoods_list) > 0:
			fragments_likelihoods_df = np.array(fragments_likelihoods_list)
			log_fragments_likelihoods_df = np.log(fragments_likelihoods_df)
			sum_log_fragments_likelihoods_df = log_fragments_likelihoods_df.sum(axis = 0)
			sum_log_fragments_likelihoods_df[np.isinf(sum_log_fragments_likelihoods_df)] = None


	return sum_log_fragments_likelihoods_df

def calculate_phred(joint_probabilities):
	libphred = np.ctypeslib.load_library('libphred.so','/groups/nshomron/tomr/projects/cffdna/hoobari/src/') #TODO - fix to relative path
	libphred.calculatePhred.restype = ctypes.c_longdouble
	cdubs = joint_probabilities.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble))

	return libphred.calculatePhred(cdubs)

def simple_phred_calculation(joint_probabilities, predicted_genotype):
	joint_probabilities_normalized = joint_probabilities - np.min(joint_probabilities[np.isfinite(joint_probabilities)])
	exp_joint_probabilities = np.exp(joint_probabilities_normalized)
	exp_joint_probabilities[np.isnan(exp_joint_probabilities)] = 0
	sum_exp_joint_probabilities = exp_joint_probabilities.sum()
	posteriors = np.asarray((exp_joint_probabilities / sum_exp_joint_probabilities), dtype = np.float)
	#print((int(exp_joint_probabilities[predicted_genotype])) / int(sum_exp_joint_probabilities))
	max_posterior = int(exp_joint_probabilities[predicted_genotype]) / int(sum_exp_joint_probabilities)


	# if predicted_genotype == 0:
	# 	phred = -10 * np.log10(1-posteriors[0])
	# elif predicted_genotype in (1, 2):
	# 	phred = -10 * np.log10(posteriors[0])

	phred = -10 * np.log10(posteriors[0])

	if np.isinf(phred):
		phred = 99999

	return phred

def calculate_posteriors(var_priors, var_likelihoods):
	# Convert to numeric values just in case
	#var_priors, var_likelihoods = pd.to_numeric(var_priors), pd.to_numeric(var_likelihoods)
	var_priors, var_likelihoods = np.array(var_priors, dtype = np.float), np.array(var_likelihoods, dtype = np.float)
	print(var_priors, var_likelihoods)
	# sum priors and likelihoods
	# if there are no priors (for instance if parental genotypes at positions are missing),
	# take only likelihoods
	joint_probabilities = np.add(var_priors, var_likelihoods)
	probabilities_source = 'joint'
	if not len(joint_probabilities[~np.isnan(joint_probabilities)]) > 0: # if there are only nans in the joint #0 IS FALSE IN ANY!!!!!! USE LEN()
		if len(var_likelihoods[~np.isnan(var_likelihoods)]) > 0:
			joint_probabilities = var_likelihoods
			probabilities_source = 'likelihoods'
		# elif any(var_priors[~np.isnan(var_priors)]):
		# 	joint_probabilities = var_priors
		# 	probabilities_source = 'priors'
		else:
			joint_probabilities = prediction = phred = None
			probabilities_source = 'none'
			
	if joint_probabilities is not None:
		joint_probabilities[np.isnan(joint_probabilities)] = -999
		prediction = joint_probabilities.argmax()
		# phred = calculate_phred(joint_probabilities)
		# phred = simple_phred_calculation(joint_probabilities, prediction)
		phred = 300

	return (joint_probabilities, prediction, phred, probabilities_source)
