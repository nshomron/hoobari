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

# project's
import parse_gt
from stderr import *
import vcfuid

# Note that all database operations should be in their respective module
# This could drastically simplify operations like multiple database handling (if we keep it)
from db import Variants

def get_fetal_and_shared_lengths(db_path):

    con = Variants(db_path, probe=False)

    shared_df = con.fetalLengthDist()
    fetal_df = con.sharedLengthDist()

    return (shared_df, fetal_df)

def create_length_distributions(db_path, cores = False, db_prefix = False):
	'''
	variant name - form of chr:position_ref/alt
	fetal - choose either the fetal allele is suppose to be the ref (if mother = 1/1 and father = 0/0)
	or the alt (if mother = 0/0 and father = 1/1)
	'''

	if cores:
		pool = Pool(int(cores))
	else:
		if cpu_count() > 1:
			pool = Pool(cpu_count() - 1)
		else:
			pool = Pool(1)

	db_path = os.path.abspath(db_path)
	if db_prefix:
		db_file_regex = re.compile(db_prefix + r'.*\.db')
		db_files_loc = os.path.dirname(db_path)
		db_files = [os.path.join(db_files_loc, file) for file in os.listdir(db_files_loc) if re.match(db_file_regex, file)]
	else:
		db_files = [db_path]

	pooled_results = pool.map(get_fetal_and_shared_lengths, db_files)
	pool.close()
	pool.join()

	shared_lengths = pooled_results[0][0]
	fetal_lengths = pooled_results[0][1]
	for tup in pooled_results[1:]:
		shared_lengths = shared_lengths.merge(tup[0], how='outer')
		fetal_lengths = fetal_lengths.merge(tup[1], how='outer')

    #TODO: Can I disable sorting?
	shared_lengths = shared_lengths.groupby(by='length').sum()
	fetal_lengths = fetal_lengths.groupby(by='length').sum()

	return (shared_lengths, fetal_lengths)

def generate_length_distributions_plot(shared_lengths, fetal_lengths, file_name):
	#show length distributions plot
	shared_lengths[shared_lengths.index < 501].plot()
	fetal_lengths[fetal_lengths.index < 501].plot()
	plt.savefig(file_name)

def create_fetal_fraction_per_length_df(shared_lengths, fetal_lengths, fetal_fraction, window = 3, max_len = 500):
	'''
	make a dictionary that shows the fetal fraction at each fragment length
	'''

	# generate fetal fraction per fragment length (window) table
	bins = range(0, max_len, window)

	shared_lengths_list = []
	for i,c in shared_lengths.iterrows():
		for j in range(int(c)):
			shared_lengths_list.append(int(i))

	fetal_lengths_list = []
	for i,c in fetal_lengths.iterrows():
		for j in range(int(c)):
			fetal_lengths_list.append(int(i))

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

	fetal_fraction_per_length_df = pd.Series(index = range(0, bins[-1] + 1, 1))

	binned_list_indices = [0] + list(nprepeat(range(len(fetal_binned)), window))

	for i in range(len(fetal_fraction_per_length_df)):
		idx_in_binned = binned_list_indices[i]
		fetal = fetal_binned[idx_in_binned][0]
		shared = shared_binned[idx_in_binned][0]
		if fetal > 10 and shared > 10:
			ff = (2 * fetal) / (shared + fetal)
			if ff > 1:
				if i > 0:
					fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
				else:
					fetal_fraction_per_length_df[i] = fetal_fraction
			else:
				fetal_fraction_per_length_df[i] = ff
		else:
			if i > 0:
				fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
			else:
				fetal_fraction_per_length_df[i] = fetal_fraction

	printverbose(fetal_fraction_per_length_df)

	return fetal_fraction_per_length_df

def calculate_total_fetal_fraction(shared_lengths, fetal_lengths):

	n_shared = int(shared_lengths.sum())
	n_fetal = int(fetal_lengths.sum())
	total_fetal_fraction = (2 * n_fetal) / (n_shared + n_fetal)
	printerr('total fetal fraction:', total_fetal_fraction)
	return total_fetal_fraction

def calculate_err_rate():
	err = 0.003
	return err

def run_full_preprocessing(db_path, cores = False, db_prefix = False, window = 3, max_len = 500, plot = False):

	printerr('pre-processing', 'creating length distributions')
	shared_lengths, fetal_lengths = create_length_distributions(db_path, cores = cores, db_prefix = db_prefix)
	printerr('pre-processing', 'calculating error rate')
	err_rate = calculate_err_rate()
	printerr('pre-processing', 'calculating total fetal fraction')
	total_fetal_fraction = calculate_total_fetal_fraction(shared_lengths, fetal_lengths)
	printerr('pre-processing', 'calculating fetal fraction per read template length')
	fetal_fractions_df = create_fetal_fraction_per_length_df(shared_lengths, fetal_lengths, total_fetal_fraction, window = window, max_len = max_len)
	if plot:
		plot_file_name = 'length_distributions.png'
		printerr('pre-processing', 'saving length distributions plot as ')
		generate_length_distributions_plot(shared_lengths, fetal_lengths, plot_file_name)

	return (err_rate, total_fetal_fraction, fetal_fractions_df)
