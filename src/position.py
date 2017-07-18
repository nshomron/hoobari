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
from decimal import Decimal, getcontext
getcontext().prec = 10000
# project's
import parse_gt
from stderr import printerr
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
	else: # get maf or af or ldaf
		maf = retrieve_maf_from_ensembl(variant_name)
		if maf[1] is not None:
			p_maternal_alt = p_paternal_alt = maf[1]
			priors_source = maf[0]
		else:
			priors = [None, None, None]
			priors_source = '.'

	priors = [	(1-p_maternal_alt)*(1-p_paternal_alt),
				p_maternal_alt*(1-p_paternal_alt) + (1-p_maternal_alt)*p_paternal_alt,
				p_maternal_alt*p_paternal_alt]			
	
	for i in range(len(priors)):
		if priors[i] == 0:
			priors[i] = None
		else:
			priors[i] = np.log(priors[i])


	return (priors, priors_source)

def calculate_fragment_i(fetal_genotypes, frag_genotype, ref, alt, f, err_rate):
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

	if frag_genotype == alt:
		frag_i_likelihoods = [0.5*(1-f), 0.5, 0.5*(1+f)]
		return frag_i_likelihoods
	elif frag_genotype == ref:
		frag_i_likelihoods = [0.5*(1+f), 0.5, 0.5*(1-f)]
		return frag_i_likelihoods

def calculate_likelihoods(
	variant,
	tmp_dir,
	total_fetal_fraction,
	fetal_fractions_df,
	err_rate,
	lengths = False,
	origin = False,
	**kwargs):
	
	'''
	calculate likelihoods for each of the 3 possible fetal genotype, based on the prbability that 
	fragment_i at the position, will show a certain allele (ref or alt), given other factors of 
	the model, such as the maternal genotype, fragment length and the fetal genotype (which is unknown,
	so we check for all possibilities - 1/1, 0/1 and 0/0)
	'''

	if origin:
		lengths = True
	
	chrom, pos, ref, alt = vcfuid.uid_to_rec(variant)
	chrom_pos = str(chrom) + ':' + str(pos)

	snp_json_path = os.path.join(tmp_dir, 'jsons', chrom + '_snps', chrom_pos + '.json')
	pos_data = pd.DataFrame(json_load(snp_json_path))

	if pos_data is not None:

		fragments_likelihoods_list = []
		for row in pos_data.itertuples():
			frag_genotype = row[1]

			# get fetal fraction depending on factors
			if origin:
				frag_qname = row[3]
				if frag_qname in known_fetal_qnames_dic:
					ff = 0.7
				else: 
					try:
						ff = fetal_fractions_df[int(row[2])]
					except:
						ff = total_fetal_fraction
			elif lengths:
				try:
					ff = fetal_fractions_df[int(row[2])]
				except:
					ff = total_fetal_fraction
			else:
				ff = total_fetal_fraction


			frag_i_likelihood_list = calculate_fragment_i(fetal_genotypes, frag_genotype, ref, alt, ff, err_rate)
			if frag_i_likelihood_list is not None:
				fragments_likelihoods_list.append(frag_i_likelihood_list)
		
		fragments_likelihoods_df = np.array(fragments_likelihoods_list)
		log_fragments_likelihoods_df = np.log(fragments_likelihoods_df)
		sum_log_fragments_likelihoods_df = log_fragments_likelihoods_df.sum(axis = 0)
	
	return sum_log_fragments_likelihoods_df

def calculate_posteriors(var_priors, var_likelihoods):
	
	# sum priors and likelihoods
	# if there are no priors (for instance if parental genotypes at positions are missing),
	# take only likelihoods
	if var_priors:
		joint_probabilities = np.sum((var_priors, var_likelihoods), axis = 0)
	else:
		joint_probabilities = var_likelihoods

	# closest to 0 is the prediction (all three fetal genotypes have same number of fragments)
	prediction = np.idxmax(joint_probabilities[~np.isnan(joint_probabilities)])


	# -------- maybe this as new function? (make_phred or something) -----------
	# calculate probabilities for the output vcf and for plotting success rates per maximum posterior threshold
	joint_probabilities_c = joint_probabilities - np.min(joint_probabilities[~np.isnan(joint_probabilities)])
	exp_joint_probabilities = np.exp(joint_probabilities_c)
	if any(exp_joint_probabilities[np.isinf(exp_joint_probabilities)]):
		use_decimal = True
		exp_joint_probabilities_list = []
		for d in joint_probabilities_c:
			if np.isnan(d):
				exp_joint_probabilities_list.append(Decimal(0))
			else:
				exp_joint_probabilities_list.append(Decimal(d).exp())
		exp_joint_probabilities = np.array(exp_joint_probabilities_list)
	else:		
		use_decimal = False
		exp_joint_probabilities[np.isnan(exp_joint_probabilities)] = 0
	

	sum_exp_joint_probabilities = exp_joint_probabilities.sum()
	posteriors = np.divide(exp_joint_probabilities, sum_exp_joint_probabilities)
	posteriors = posteriors.astype(np.float128)

	# 2017-07-17: note that here there should be another condition:
	# if the max posterior shows fetal gt of 0, then phred = -10*log(P(gt is 1) + P(gt is 2)) = -10log(1 - P(gt is 0))
	# if the max posterior shows fetal gt of either 1 or 2, then phred = -10*log(1 - (P(gt is 1) + P(gt is 2))) = -10log(P(gt is 0))

	if use_decimal: 
		phred = float(-10*(1 - posteriors.max()).log10())
	else:
		phred = float(-10*(np.log10(1 - np.max(posteriors))))
	
	# ------------------------------------------------------------------------------

	return (joint_probabilities, prediction, phred)