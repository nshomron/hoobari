# --------- import modules ------------
import re
import os
import sys
import subprocess
import requests
import vcf
import numpy as np
import pandas as pd
import seaborn as sns
import json
import argparse
from multiprocessing import Pool, cpu_count
import parse_gt
from stderr import printerr

# --------- functions ----------

def json_load(path):
	with open(path, 'r') as f:
		json_object = json.load(f)
	return json_object

def get_qnames_and_alleles(variant_name):
	pos_ref_frags_dic = {}
	pos_alt_frags_dic = {}
	
	var_split = re.split(r':|_|/', variant_name)
	# print(var_split)
	
	chrom = var_split[0]
	chrom_pos = var_split[0] + ':' + str(var_split[1])
	ref = var_split[2]
	alt = var_split[3]

	try:
		snp_json_path = os.path.join(args.tmp_dir, 'jsons', chrom + '_snps', chrom_pos + '.json')
		pos_data = pd.DataFrame(json_load(snp_json_path))
		# print(pos_data)
		ref_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == ref].iloc[:,1:3]
		for row in ref_lengths_and_qnames_at_position_df.itertuples():
			pos_ref_frags_dic[row[2]] = row[1]
		
		alt_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == alt].iloc[:,1:3]
		for row in alt_lengths_and_qnames_at_position_df.itertuples():
			pos_alt_frags_dic[row[2]] = row[1]
	except:
		pass

	return pos_ref_frags_dic, pos_alt_frags_dic

def create_length_distributions(variant_list):
	'''
	variant name - form of chr:position_ref/alt
	fetal - choose either the fetal allele is suppose to be the ref (if mother = 1/1 and father = 0/0)
	or the alt (if mother = 0/0 and father = 1/1)
	'''
	
	from multiprocessing import Pool, cpu_count
	pool = Pool(cpu_count() - 1)
	pooled_results = pool.map(get_qnames_and_alleles, variant_list)

	return pooled_results

def create_fetal_fraction_per_length_df(fetal_lengths_list, shared_lengths_list, window = 3, max = 500):
	bins = range(0,max,window)
	
	fetal_pd_cut = pd.cut(fetal_lengths_list, bins, include_lowest = True)
	shared_pd_cut = pd.cut(shared_lengths_list, bins, include_lowest = True)

	
	fetal_binned = pd.value_counts(fetal_pd_cut, sort = False).to_frame().values.tolist()
	shared_binned = pd.value_counts(shared_pd_cut, sort = False).to_frame().values.tolist()
	
	
	fetal_fraction_per_length_df = pd.Series(index = range(0, bins[-1] + 1, 1))

	binned_list_indices = [0] + list(np.repeat(range(len(fetal_binned)), window))

	for i in range(len(fetal_fraction_per_length_df)):
		idx_in_binned = binned_list_indices[i]
		fetal = fetal_binned[idx_in_binned][0]
		shared = shared_binned[idx_in_binned][0]
		if fetal > 5 and shared > 5:
			fetal_fraction_per_length_df[i] = (2 * fetal) / (shared + fetal)
		else:
			fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
		#stderr(fetal_fraction_per_length_df[i])
	return fetal_fraction_per_length_df


# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_id", help = '')
parser.add_argument("-parents", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-tmp", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-no_pkls", "--no_pkls", action = 'store_true', help = "don't make or use pickle files, usually after loading large objects")
parser.add_argument("-v", "--vcf_output", default = 'stdout', help = 'path for vcf output')

args = parser.parse_args()

###print the command with all arguments


# --------- constants ----------
sample_id = args.sample_id
project_dir = os.path.join(os.path.expanduser('~/nshomron_tomr/projects/cffdna/simulations/lo/lo_from_tmp'), sample_id, 'wgs')
parents_vcf_file = os.path.join(project_dir, sample_id + '_parents_wgs.vcf')
fb_by_chrom_dir = os.path.join(project_dir, 'fb_by_chrom')
error_rate_cfdna_vcf = os.path.join(project_dir, sample_id + '_cfdna_error_positions.vcf')

# --------- main ----------
# parse vcf files
# parents
parents_gt = parse_gt.parse_split(project_dir, '.*parents.*chr.*vcf')
parents_gt = parents_gt.drop([var for var in parents_gt.index.values if var.startswith('chrY')], axis=0)
parents_gt = parents_gt.dropna()
# maternal plasma
cfdna_gt = parse_gt.parse_split(fb_by_chrom_dir, '.*cfdna.*chr.*vcf')
cfdna_gt = cfdna_gt.drop([var for var in cfdna_gt.index.values if var.startswith('chrY')], axis=0)
cfdna_gt = cfdna_gt.dropna()
# intersect to get positions of interest. this shouldn't drop to many positions, since variant calling was done by same coordinates.
positions_to_predict = set(parents_gt.index.values).intersection(cfdna_gt.index.values)
positions_to_predict = set(parents_gt.query('(M' + sample_id + 'W == 1)').index.values).intersection(cfdna_gt.index.values)

# calculations of fragment length distributions, and fetal fraction
stderr('creating lengths distributions at mother_ref_father_alt...')
mother_ref_father_alt = parents_gt.query('(M' + sample_id + 'W == 0) & (H' + sample_id + 'W == 2)').index.values
mother_ref_father_alt_dics_list = create_length_distributions(mother_ref_father_alt)
stderr('creating lengths distributions at mother_alt_father_ref...')
mother_alt_father_ref = parents_gt.query('(M' + sample_id + 'W == 2) & (H' + sample_id + 'W == 0)').index.values
mother_alt_father_ref_dics_list = create_length_distributions(mother_alt_father_ref)
stderr('creating lengths distributions dictionaries...')
shared_fragments_dic = {}
fetal_fragments_dic = {}
for tup in mother_ref_father_alt_dics_list:
	shared_fragments_dic.update(tup[0])
	fetal_fragments_dic.update(tup[1])
del mother_ref_father_alt_dics_list
for tup in mother_alt_father_ref_dics_list:
	shared_fragments_dic.update(tup[1])
	fetal_fragments_dic.update(tup[0])
del mother_alt_father_ref_dics_list

total_fetal_fraction = (2 * len(fetal_fragments_dic)) / (len(shared_fragments_dic) + len(fetal_fragments_dic))
stderr('Total fetal fraction: ', total_fetal_fraction)

shared_fragment_lengths_list = list(shared_fragments_dic.values())
fetal_fragment_lengths_list = list(fetal_fragments_dic.values())

# show length distributions plot
# shared_plot_data = [i for i in shared_fragment_lengths_list if i < 500]
# fetal_plot_data = [i for i in fetal_fragment_lengths_list if i < 500]
# sns.kdeplot(np.array(shared_plot_data), bw = 0.15)
# sns.set_style('whitegrid')
# sns.kdeplot(np.array(fetal_plot_data), bw = 0.15)
# sns.plt.show()

# creating a pool of known fetal fragments to add as data to the model
stderr('creating a list of known fetal fragments...')
stderr('extracting fetal fragments from mother_ref_father_het positions...')
mother_ref_father_het = parents_gt.query('(M' + sample_id + 'W == 0) & (H' + sample_id + 'W == 1)').index.values
mother_ref_father_het_dics_list = create_length_distributions(mother_ref_father_het)
for tup in mother_ref_father_het_dics_list:
	fetal_fragments_dic.update(tup[1])
del mother_ref_father_het_dics_list
stderr('extracting fetal fragments from mother_alt_father_het positions...')
mother_alt_father_het = parents_gt.query('(M' + sample_id + 'W == 2) & (H' + sample_id + 'W == 1)').index.values
mother_alt_father_het_dics_list = create_length_distributions(mother_alt_father_het, fb_by_chrom_dir, fetal = 'ref')
for tup in mother_alt_father_het_dics_list:
	fetal_fragments_dic.update(tup[0])
del mother_alt_father_het_dics_list

# make a dictionary that shows the fetal fraction at each fragment length
stderr('making a dictionary that shows the fetal fraction at each fragment length:')
ff_per_length_df = create_fetal_fraction_per_length_df(fetal_fragment_lengths_list, shared_fragment_lengths_list)