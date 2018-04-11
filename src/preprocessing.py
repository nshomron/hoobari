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
import db


# --------- functions ----------

def get_fetal_and_shared_lengths(db_path, qnames = False):
	'''
	input - path to the database from hoobari's patch
	output - a tuple with two dictionaries, one contains the counts of fetal fragments at different lengths,
	and the other is similar, but for fragments which aren't necessarily fetal ("shared")
	'''
	con = db.Variants(db_path, probe=False)
	fetal_lengths = con.getFetalLengths()
	shared_lengths = con.getSharedLengths()

	if qnames:
		fetal_qnames, shared_qnames = con.getFetalSharedQnames()
		return (shared_lengths, fetal_lengths, shared_qnames, fetal_qnames)
	else:
		return (shared_lengths, fetal_lengths)


def create_length_distributions(db_path, cores, db_prefix = False, qnames = False):
	'''
	input: db_path - path to the database from hoobari's patch; cores - number of cores to use for
	multiprocessing; db_prefix - if hoobari's patch was ran for many different regions, it creates
	many DB's with the same prefix. since all those databases are required in this step, the prefix
	is also required.
	output: a tuple with two pandas dataframes that contain fragment length distributions which were
	calculated using *all* the databases
	'''

	# if the number of cores was specified - use it. if not, and there's more than 1 core, use all cores
	# except for one.
	# if cores:
	# 	pool = Pool(int(cores))
	# else:
	# 	if cpu_count() > 1:
	# 		pool = Pool(cpu_count() - 1)
	# 	else:
	# 		pool = Pool(1)
	pool = Pool(int(cores))

	# if db_prefix was given as an argument, work for all databases with that prefix.
	db_path = os.path.abspath(db_path)
	if db_prefix:
		db_file_regex = re.compile(db_prefix + r'.*\.db')
		db_files_loc = os.path.dirname(db_path)
		db_files = [os.path.join(db_files_loc, file) for file in os.listdir(db_files_loc) if re.match(db_file_regex, file)]
	else:
		db_files = [db_path]

	# run the function get_fetal_and_shared_lengths for each path in db_files
	# pooled_results = pool.map(get_fetal_and_shared_lengths, db_files)
	get_qnames_and_alleles_with_args = partial(get_fetal_and_shared_lengths, qnames = qnames)
	# pooled_results = pool.map(get_qnames_and_alleles_with_args, db_files) #TODO: use pool.imap_unordered(func, iterable[, chunksize])
	# pool.close()
	# pool.join()

	# create two lists, one with all the shared fragments results, and one for the fetal fragments results
	
	con = db.Variants(db_files[0], probe=False)
	shared_lengths = con.getSharedLengths()
	fetal_lengths = con.getFetalLengths()
	if qnames:
		fetal_qnames, shared_qnames = con.getFetalSharedQnames()

	for tup in pool.imap_unordered(get_qnames_and_alleles_with_args, db_files[1:]):
		shared_lengths = shared_lengths.add(tup[0], fill_value=0)
		fetal_lengths = fetal_lengths.add(tup[1], fill_value=0)
		if qnames:
			shared_qnames_set.update(tup[2])
			fetal_qnames_set.update(tup[3])

	# for db_path in db_files[1:]:
	# 	con = db.Variants(db_path, probe=False)
	# 	shared_lengths = shared_lengths.add(con.getSharedLengths(), fill_value=0)
	# 	fetal_lengths = fetal_lengths.add(con.getFetalLengths(), fill_value=0)
	# 	if qnames:
	# 		tup = con.getFetalSharedQnames()
	# 		shared_qnames_set.update(tup[1])
	# 		fetal_qnames_set.update(tup[0])

	if qnames:
		with open('shared_qnames_list.txt', 'w') as f:
			for q in shared_qnames_set:
				print(q, file = f)
		with open('fetal_qnames_list.txt', 'w') as f:
			for q in fetal_qnames_set:
				print(q, file = f)

	return(shared_lengths, fetal_lengths)

def estimate_length_distribution(total_fetal_fraction):
	normalized_fetal_fractions_dist = pd.Series([1.195695703013115, 1.195695703013115, 1.195695703013115, 1.195695703013115, 2.1292505317701584, 2.1292505317701584, 2.1292505317701584, 2.266846530994138, 2.266846530994138, 2.266846530994138, 1.500350157313987, 1.500350157313987, 1.500350157313987, 1.4017112196359047, 1.4017112196359047, 1.4017112196359047, 2.1773657468759486, 2.1773657468759486, 2.1773657468759486, 1.2377350652159864, 1.2377350652159864, 1.2377350652159864, 1.5465401566595522, 1.5465401566595522, 1.5465401566595522, 2.4365123449703865, 2.4365123449703865, 2.4365123449703865, 1.4853524240997755, 1.4853524240997755, 1.4853524240997755, 1.3845453630979483, 1.3845453630979483, 1.3845453630979483, 1.502413694347957, 1.502413694347957, 1.502413694347957, 1.4704651566050584, 1.4704651566050584, 1.4704651566050584, 1.3658814263443844, 1.3658814263443844, 1.3658814263443844, 1.3763656011350447, 1.3763656011350447, 1.3763656011350447, 1.4051395839047558, 1.4051395839047558, 1.4051395839047558, 1.418233691511955, 1.418233691511955, 1.418233691511955, 1.3090329283666646, 1.3090329283666646, 1.3090329283666646, 1.431968780624545, 1.431968780624545, 1.431968780624545, 1.5523974278336532, 1.5523974278336532, 1.5523974278336532, 1.5055667413383695, 1.5055667413383695, 1.5055667413383695, 1.5402599292330517, 1.5402599292330517, 1.5402599292330517, 1.6017058738173935, 1.6017058738173935, 1.6017058738173935, 1.8423934249381364, 1.8423934249381364, 1.8423934249381364, 1.5959188109726787, 1.5959188109726787, 1.5959188109726787, 1.6819447781494345, 1.6819447781494345, 1.6819447781494345, 1.838203977414241, 1.838203977414241, 1.838203977414241, 1.7224392310333196, 1.7224392310333196, 1.7224392310333196, 1.7289746387048759, 1.7289746387048759, 1.7289746387048759, 1.9095170280669445, 1.9095170280669445, 1.9095170280669445, 2.072064953768226, 2.072064953768226, 2.072064953768226, 2.031656551788986, 2.031656551788986, 2.031656551788986, 2.0585066571357045, 2.0585066571357045, 2.0585066571357045, 2.3809272197788394, 2.3809272197788394, 2.3809272197788394, 2.1077576575995667, 2.1077576575995667, 2.1077576575995667, 2.191545927606252, 2.191545927606252, 2.191545927606252, 2.2134804031902116, 2.2134804031902116, 2.2134804031902116, 2.265610562618742, 2.265610562618742, 2.265610562618742, 2.2168216834122103, 2.2168216834122103, 2.2168216834122103, 2.218832271923642, 2.218832271923642, 2.218832271923642, 2.2368565308945136, 2.2368565308945136, 2.2368565308945136, 2.25614804726495, 2.25614804726495, 2.25614804726495, 2.254777147721188, 2.254777147721188, 2.254777147721188, 2.174109728454312, 2.174109728454312, 2.174109728454312, 2.0921539080118357, 2.0921539080118357, 2.0921539080118357, 1.9909559873678913, 1.9909559873678913, 1.9909559873678913, 1.8062678281374691, 1.8062678281374691, 1.8062678281374691, 1.7368734604292917, 1.7368734604292917, 1.7368734604292917, 1.6218278923731848, 1.6218278923731848, 1.6218278923731848, 1.4322992178660081, 1.4322992178660081, 1.4322992178660081, 1.3352647823898771, 1.3352647823898771, 1.3352647823898771, 1.2686590417033243, 1.2686590417033243, 1.2686590417033243, 1.0927236004122058, 1.0927236004122058, 1.0927236004122058, 0.9803002372280191, 0.9803002372280191, 0.9803002372280191, 0.8427104331050591, 0.8427104331050591, 0.8427104331050591, 0.7381302643241816, 0.7381302643241816, 0.7381302643241816, 0.7075162035381661, 0.7075162035381661, 0.7075162035381661, 0.676098160742928, 0.676098160742928, 0.676098160742928, 0.6278752302805902, 0.6278752302805902, 0.6278752302805902, 0.6231931763990364, 0.6231931763990364, 0.6231931763990364, 0.63257035977172, 0.63257035977172, 0.63257035977172, 0.6139473479557893, 0.6139473479557893, 0.6139473479557893, 0.5857105932383947, 0.5857105932383947, 0.5857105932383947, 0.5881352460327547, 0.5881352460327547, 0.5881352460327547, 0.594082727260782, 0.594082727260782, 0.594082727260782, 0.586817647502322, 0.586817647502322, 0.586817647502322, 0.5698097918562272, 0.5698097918562272, 0.5698097918562272, 0.5806322325263193, 0.5806322325263193, 0.5806322325263193, 0.5904635937592888, 0.5904635937592888, 0.5904635937592888, 0.5840741907923023, 0.5840741907923023, 0.5840741907923023, 0.5840255259523812, 0.5840255259523812, 0.5840255259523812, 0.5922311565838844, 0.5922311565838844, 0.5922311565838844, 0.6064678225478712, 0.6064678225478712, 0.6064678225478712, 0.6136005641432725, 0.6136005641432725, 0.6136005641432725, 0.6244629891099108, 0.6244629891099108, 0.6244629891099108, 0.6376257216882676, 0.6376257216882676, 0.6376257216882676, 0.6542991471721685, 0.6542991471721685, 0.6542991471721685, 0.7033573442096076, 0.7033573442096076, 0.7033573442096076, 0.703911126676141, 0.703911126676141, 0.703911126676141, 0.7587414061130653, 0.7587414061130653, 0.7587414061130653, 0.8740076538936232, 0.8740076538936232, 0.8740076538936232, 0.884150425181996, 0.884150425181996, 0.884150425181996, 0.9426520776310289, 0.9426520776310289, 0.9426520776310289, 1.1282294239906596, 1.1282294239906596, 1.1282294239906596, 1.2153840648878125, 1.2153840648878125, 1.2153840648878125, 1.2079720982674649, 1.2079720982674649, 1.2079720982674649, 1.4132332099959584, 1.4132332099959584, 1.4132332099959584, 1.6230888667079653, 1.6230888667079653, 1.6230888667079653, 1.6410133756480851, 1.6410133756480851, 1.6410133756480851, 1.5881958805960414, 1.5881958805960414, 1.5881958805960414, 1.804646175017287, 1.804646175017287, 1.804646175017287, 1.8663031895199826, 1.8663031895199826, 1.8663031895199826, 1.722164272055739, 1.722164272055739, 1.722164272055739, 1.647652996618907, 1.647652996618907, 1.647652996618907, 1.6469506154771867, 1.6469506154771867, 1.6469506154771867, 1.5591307817064606, 1.5591307817064606, 1.5591307817064606, 1.4719108420425855, 1.4719108420425855, 1.4719108420425855, 1.4297892857851822, 1.4297892857851822, 1.4297892857851822, 1.2920953050015427, 1.2920953050015427, 1.2920953050015427, 1.2205711342919268, 1.2205711342919268, 1.2205711342919268, 1.1685208246563323, 1.1685208246563323, 1.1685208246563323, 1.049551130346437, 1.049551130346437, 1.049551130346437, 0.9610278990446985, 0.9610278990446985, 0.9610278990446985, 0.8977333939331735, 0.8977333939331735, 0.8977333939331735, 0.8334190606426094, 0.8334190606426094, 0.8334190606426094, 0.7648045083513764, 0.7648045083513764, 0.7648045083513764, 0.6988298101538835, 0.6988298101538835, 0.6988298101538835, 0.6585666878492391, 0.6585666878492391, 0.6585666878492391, 0.6111684996680828, 0.6111684996680828, 0.6111684996680828, 0.5665580244134178, 0.5665580244134178, 0.5665580244134178, 0.5148781211828849, 0.5148781211828849, 0.5148781211828849, 0.48618640994323403, 0.48618640994323403, 0.48618640994323403, 0.46033587532991155, 0.46033587532991155, 0.46033587532991155, 0.42457852975307236, 0.42457852975307236, 0.42457852975307236, 0.3910311549657173, 0.3910311549657173, 0.3910311549657173, 0.36425156818383675, 0.36425156818383675, 0.36425156818383675, 0.3470777453124518, 0.3470777453124518, 0.3470777453124518, 0.3297880637926003, 0.3297880637926003, 0.3297880637926003, 0.3107121988327585, 0.3107121988327585, 0.3107121988327585, 0.29652484531889745, 0.29652484531889745, 0.29652484531889745, 0.2811141469156768, 0.2811141469156768, 0.2811141469156768, 0.2667474285969589, 0.2667474285969589, 0.2667474285969589, 0.2617349980954065, 0.2617349980954065, 0.2617349980954065, 0.24441903266922857, 0.24441903266922857, 0.24441903266922857, 0.23718219453784634, 0.23718219453784634, 0.23718219453784634, 0.23607152333768355, 0.23607152333768355, 0.23607152333768355, 0.23895766208345, 0.23895766208345, 0.23895766208345, 0.23069303775657976, 0.23069303775657976, 0.23069303775657976, 0.22778154338144438, 0.22778154338144438, 0.22778154338144438, 0.22227276488421002, 0.22227276488421002, 0.22227276488421002, 0.2259470817454212, 0.2259470817454212, 0.2259470817454212, 0.2286492302134698, 0.2286492302134698, 0.2286492302134698, 0.23503432221840026, 0.23503432221840026, 0.23503432221840026, 0.23077484051921943, 0.23077484051921943, 0.23077484051921943, 0.2422018108014714, 0.2422018108014714, 0.2422018108014714, 0.2406218002553601, 0.2406218002553601, 0.2406218002553601, 0.26980269067502727, 0.26980269067502727, 0.26980269067502727, 0.2723817890491493, 0.2723817890491493, 0.2723817890491493, 0.29331623708078347, 0.29331623708078347, 0.29331623708078347, 0.33219796463591716, 0.33219796463591716, 0.33219796463591716, 0.3292209097084398, 0.3292209097084398, 0.3292209097084398, 0.39068358636774314, 0.39068358636774314, 0.39068358636774314, 0.45534117696074344, 0.45534117696074344, 0.45534117696074344, 0.49527403960628974, 0.49527403960628974, 0.49527403960628974, 0.5293983772269767, 0.5293983772269767, 0.5293983772269767, 0.5805258692936526, 0.5805258692936526, 0.5805258692936526, 0.5972556785462638, 0.5972556785462638, 0.5972556785462638, 0.6297278740840514, 0.6297278740840514, 0.6297278740840514, 0.6637066874967081, 0.6637066874967081, 0.6637066874967081, 0.7017438932271006, 0.7017438932271006, 0.7017438932271006, 0.6973996643075412, 0.6973996643075412, 0.6973996643075412, 0.6460060403539388, 0.6460060403539388, 0.6460060403539388, 0.6474303229677442, 0.6474303229677442, 0.6474303229677442, 0.6356648095891243, 0.6356648095891243, 0.6356648095891243, 0.5957781895580423, 0.5957781895580423, 0.5957781895580423, 0.5625754979314425, 0.5625754979314425, 0.5625754979314425, 0.5474479316706278, 0.5474479316706278, 0.5474479316706278, 0.5020326892714471, 0.5020326892714471, 0.5020326892714471, 0.4812812131075173, 0.4812812131075173, 0.4812812131075173, 0.4625050816019519, 0.4625050816019519, 0.4625050816019519, 0.4551534498168077, 0.4551534498168077, 0.4551534498168077, 0.41882891956888196, 0.41882891956888196, 0.41882891956888196, 0.40052291725759387, 0.40052291725759387, 0.40052291725759387, 0.3989206002929431, 0.3989206002929431, 0.3989206002929431, 0.3844624865835286, 0.3844624865835286, 0.3844624865835286, 0.3379402117997727, 0.3379402117997727, 0.3379402117997727])
	fetal_fractions_df = total_fetal_fraction * normalized_fetal_fractions_dist
	return(fetal_fractions_df)

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

	# n_shared = int(shared_lengths[shared_lengths.index < 501].sum())
	# n_fetal = int(fetal_lengths[fetal_lengths.index < 501].sum())
	n_shared = int(shared_lengths.sum())
	n_fetal = int(fetal_lengths.sum())
	total_fetal_fraction = (2 * n_fetal) / (n_shared + n_fetal)
	
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
							total_fetal_fraction,
							calculate_empirical_ff_dist,
							cores,
							db_prefix = False,
							window = False,
							max_len = 500,
							plot = False,
							qnames = False):

	printerr('pre-processing', 'calculating error rate (a place holder is temporarily set to 0.003)')
	err_rate = calculate_err_rate()

	calculate_fetal_fraction = not total_fetal_fraction
	
	if calculate_empirical_ff_dist or calculate_fetal_fraction:
		printerr('pre-processing', 'creating length distributions')
		shared_lengths, fetal_lengths = create_length_distributions(db_path, cores, db_prefix = db_prefix, qnames = qnames)
		if plot:
			printerr('pre-processing', 'saving length distributions plot as', fetal_sample + '.length_distributions.png')
			generate_length_distributions_plot(shared_lengths, fetal_lengths, fetal_sample)
		if calculate_fetal_fraction:
			printerr('pre-processing', 'calculating total fetal fraction')
			total_fetal_fraction = calculate_total_fetal_fraction(shared_lengths, fetal_lengths)
	
	printerr('total fetal fraction:', total_fetal_fraction)

	if calculate_empirical_ff_dist:
		printerr('pre-processing', 'calculating an empirical distribution of fetal fractions per fragment length')
		fetal_fractions_df = create_fetal_fraction_per_length_df(shared_lengths,
									fetal_lengths,
									total_fetal_fraction,
									err_rate,
									window,
									max_len = max_len)
	else:
		printerr('calculating an estimated distribution of fetal fractions per fragment length based on prior knowledge')
		fetal_fractions_df = estimate_length_distribution(total_fetal_fraction)

		

	return (err_rate, total_fetal_fraction, fetal_fractions_df)
