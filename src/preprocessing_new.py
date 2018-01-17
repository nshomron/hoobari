# --------- import modules ------------
# external
import os, sys
import sqlite3
import vcf
from collections import Counter
from numpy import repeat as nprepeat
import pandas as pd
import matplotlib
matplotlib.use('Agg') # to allow saving a figure even though display is invalid
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
from functools import partial
import re
import db

# project's
import parse_gt
from stderr import *
import vcfuid



# --------- functions ----------

def get_fetal_and_shared_lengths(db_path, qnames = False):
	'''
	input - path to the database from hoobari's patch
	output - a tuple with two dictionaries, one contains the counts of fetal fragments at different lengths,
	and the other is similar, but for fragments which aren't necessarily fetal ("shared")
	'''

	con = Variants(db_path, probe=False)
	shared_lengths, fetal_lengths = con.getSharedLengths(), con.getFetalLengths()

	if qnames:
		fetal_qnames, shared_qnames = con.getFetalSharedQnames()

		return (shared_lengths, fetal_lengths, shared_qnames, fetal_qnames)
	else:
		return (shared_lengths, fetal_lengths)

def generate_length_distributions_plot(shared_lengths, fetal_lengths, fetal_sample):
	'''
	save the length distributions plot
	'''
	fetal_lengths_500 = fetal_lengths[fetal_lengths.index < 501]
	shared_lengths_500 = shared_lengths[shared_lengths.index < 501]
	shared_density = shared_lengths_500 / shared_lengths_500.sum()
	fetal_density = fetal_lengths_500 / fetal_lengths_500.sum()
	length_distributions_df = pd.concat([fetal_density.iloc[1:500], shared_density.iloc[1:500]], axis = 1)
	length_distributions_df.columns = ['fetal fragments', 'shared fragments']
	length_distributions_df.plot()
	plt.savefig(fetal_sample + '.length_distributions.png')

def calculate_total_fetal_fraction(shared_lengths, fetal_lengths):

	n_shared = int(shared_lengths.sum())
	n_fetal = int(fetal_lengths.sum())
	total_fetal_fraction = (2 * n_fetal) / (n_shared + n_fetal)
	printerr('total fetal fraction:', total_fetal_fraction)

	return total_fetal_fraction

def create_fetal_fraction_per_length_df(shared_lengths, # pandas dataframe of shared fragments length distribution
					fetal_lengths, # pandas dataframe of fetal fragments length distribution
					total_fetal_fraction, # the total calculated percent of fetal DNA within the cfDNA
					err_rate, # TODO: will be later added to the model
					window = False, # for the length distributions. for example: 10 will give
					max_len = 500):
	'''
	output:
	make a dictionary that shows the fetal fraction at each fragment length
	input:
	shared_lengths - pandas dataframe of shared fragments length distribution
	fetal_lengths - pandas dataframe of fetal fragments length distribution
	total_fetal_fraction - the total calculated percent of fetal DNA within the cfDNA
	err_rate - TODO: will be later added to the model
	window - for the length distributions. for example: window=10 will give the same fetal fraction for fragments
	with length 0-10, 11-20, ...
	max_len - maximum fragment length used. default: 500.

	'''
	# TODO: rewrite a simpler version of this function


	# define bins (aka categories)
	bins = range(0, max_len, window)

	# make lists from two length distribution dataframes
	shared_lengths_list = []
	for i,c in shared_lengths.iterrows():
		for j in range(int(c)):
			shared_lengths_list.append(int(i))

	fetal_lengths_list = []
	for i,c in fetal_lengths.iterrows():
		for j in range(int(c)):
			fetal_lengths_list.append(int(i))

	# order the list using the bins - from a continuous variable to a categorical variable
	shared_pd_cut = pd.cut(shared_lengths_list, bins, include_lowest = True)
	fetal_pd_cut = pd.cut(fetal_lengths_list, bins, include_lowest = True)
	printverbose('shared_pd_cut')
	printverbose(shared_pd_cut)
	printverbose('fetal_pd_cut')
	printverbose(fetal_pd_cut)

	fetal_binned = pd.value_counts(fetal_pd_cut, sort = False).to_frame().values.tolist()
	shared_binned = pd.value_counts(shared_pd_cut, sort = False).to_frame().values.tolist()
	printverbose('fetal_binned')
	printverbose(fetal_binned)
	printverbose('shared_binned')
	printverbose(shared_binned)

	# for each bin, calculate its fetal fraction
	# if the fetal fraction is above 1 use (1-err);
	# if there is not enough fragments (under 5) use the result of the former window
	fetal_fraction_per_length_df = pd.Series(index = range(0, bins[-1] + 1, 1))

	binned_list_indices = [0] + list(nprepeat(range(len(fetal_binned)), window))

	for i in range(len(fetal_fraction_per_length_df)):
		idx_in_binned = binned_list_indices[i]
		fetal = fetal_binned[idx_in_binned][0]
		shared = shared_binned[idx_in_binned][0]
		if fetal > 5 or shared > 5: # TODO: or? and?
			ff = (2 * fetal) / (shared + fetal)
			if ff > 1:
				fetal_fraction_per_length_df[i] = 1 - err_rate
				# if i > 0:
				# 	fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
				# else:
				# 	fetal_fraction_per_length_df[i] = fetal_fraction
			else:
				fetal_fraction_per_length_df[i] = ff
		else:
			if i > 0:
				fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
			else:
				fetal_fraction_per_length_df[i] = total_fetal_fraction



	printverbose(fetal_fraction_per_length_df)

	return fetal_fraction_per_length_df

def calculate_err_rate():
	err = 0.003
	return err

def run_full_preprocessing(	db_path,
							fetal_sample,
							cores = False,
							db_prefix = False,
							window = False,
							max_len = 500,
							plot = False,
							qnames = False):

	printerr('pre-processing', 'creating length distributions')
	shared_lengths, fetal_lengths = create_length_distributions(db_path, cores = cores, db_prefix = db_prefix, qnames = qnames)
	printerr('pre-processing', 'calculating error rate')
	err_rate = calculate_err_rate()
	printerr('pre-processing', 'calculating total fetal fraction')
	total_fetal_fraction = calculate_total_fetal_fraction(shared_lengths, fetal_lengths)
	printerr('pre-processing', 'calculating fetal fraction per read template length')
	fetal_fractions_df = create_fetal_fraction_per_length_df(shared_lengths, fetal_lengths, total_fetal_fraction, err_rate, window, max_len = max_len)
	if plot:
		printerr('pre-processing', 'saving length distributions plot as', fetal_sample + '.length_distributions.csv')
		generate_length_distributions_plot(shared_lengths, fetal_lengths, fetal_sample)

	return (err_rate, total_fetal_fraction, fetal_fractions_df)