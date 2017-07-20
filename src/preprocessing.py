# --------- import modules ------------
# external
from os import path
import vcf
from numpy import repeat as nprepeat
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

# project's
from json_commands import *
import parse_gt
from stderr import printerr
import varuid



# --------- functions ----------

def get_qnames_and_alleles(variant_name, tmp_dir):
	pos_ref_frags_dic = {}
	pos_alt_frags_dic = {}
	
	chrom, pos, ref, alt = varuid.uid_to_rec(variant_name)

	snp_json_path = path.join(tmp_dir, 'jsons', chrom + '_snps', str(pos) + '.json')
	if path.isfile(snp_json_path):
		pos_data = pd.DataFrame(json_load(snp_json_path))

		ref_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == ref].iloc[:,1:3]
		for row in ref_lengths_and_qnames_at_position_df.itertuples():
			pos_ref_frags_dic[row[2]] = row[1]
		
		alt_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == alt].iloc[:,1:3]
		for row in alt_lengths_and_qnames_at_position_df.itertuples():
			pos_alt_frags_dic[row[2]] = row[1]

	return pos_ref_frags_dic, pos_alt_frags_dic

def create_length_distributions(variant_list):
	'''
	variant name - form of chr:position_ref/alt
	fetal - choose either the fetal allele is suppose to be the ref (if mother = 1/1 and father = 0/0)
	or the alt (if mother = 0/0 and father = 1/1)
	'''	
	
	if cpu_count() > 1:
		pool = Pool(cpu_count() - 1)
	else:
		pool = Pool(1)

	pooled_results = pool.map(get_qnames_and_alleles, variant_list)

	return pooled_results

def generate_length_distributions_plot(shared_frags_length_list, fetal_frags_length_list, directory):
	#show length distributions plot
	shared_plot_data = pd.Series(shared_frags_length_list)
	fetal_plot_data = pd.Series(fetal_frags_length_list)
	shared_plot_data[shared_plot_data < 500].plot()
	fetal_plot_data[fetal_plot_data < 500].plot()
	plt.savefig(directory, 'plots', 'length_distributions.png')

def create_fetal_fraction_per_length_df(shared_fragments_dic, fetal_fragments_dic, window = 3, max = 500, plot_dir = False):
	'''
	make a dictionary that shows the fetal fraction at each fragment length
	'''

	stderr('making a dictionary that shows the fetal fraction at each fragment length:')

	shared_fragment_lengths_list = list(shared_fragments_dic.values())
	fetal_fragment_lengths_list = list(fetal_fragments_dic.values())

	# generate length distributions plot
	if plot_dir:
		generate_length_distributions_plot(shared_fragment_lengths_list, fetal_fragment_lengths_list, plot_dir)

	# generate fetal fraction per fragment length (window) table
	bins = range(0,max,window)
	
	fetal_pd_cut = pd.cut(fetal_fragment_lengths_list, bins, include_lowest = True)
	shared_pd_cut = pd.cut(shared_fragment_lengths_list, bins, include_lowest = True)

	
	fetal_binned = pd.value_counts(fetal_pd_cut, sort = False).to_frame().values.tolist()
	shared_binned = pd.value_counts(shared_pd_cut, sort = False).to_frame().values.tolist()
	
	
	fetal_fraction_per_length_df = pd.Series(index = range(0, bins[-1] + 1, 1))

	binned_list_indices = [0] + list(nprepeat(range(len(fetal_binned)), window))

	for i in range(len(fetal_fraction_per_length_df)):
		idx_in_binned = binned_list_indices[i]
		fetal = fetal_binned[idx_in_binned][0]
		shared = shared_binned[idx_in_binned][0]
		if fetal > 5 and shared > 5:
			fetal_fraction_per_length_df[i] = (2 * fetal) / (shared + fetal)
		else:
			fetal_fraction_per_length_df[i] = fetal_fraction_per_length_df[i-1]
	
	return fetal_fraction_per_length_df

def parse_parents_vcf(parents_vcf_file, fetal_sex = None):
	parents_genotypes = parse_gt.parse(parents_vcf_file)
	if fetal_sex == 'female':
		parents_genotypes = parents_genotypes.drop([var for var in parents_genotypes.index.values if var.startswith('chrY')], axis=0)
	return parents_genotypes

def create_fetal_and_shared_fragment_pools(parents_gt, maternal_sample_name, paternal_sample_name):
	'''
	calculations of fragment length distributions
	'''

	shared_fragments_dic = {}
	fetal_fragments_dic = {}

	stderr('creating lengths distributions at mother_ref_father_alt...')
	mother_ref_father_alt = parents_gt.query(	'(' + maternal_sample_name + ' == 0) & \
							(' + paternal_sample_name + ' == 2) & \
							(indel == 0)').index.values
	mother_ref_father_alt_dics_list = create_length_distributions(mother_ref_father_alt)
	for tup in mother_ref_father_alt_dics_list:
		shared_fragments_dic.update(tup[0]) # ref is shared
		fetal_fragments_dic.update(tup[1]) # alt is fetal
	del mother_ref_father_alt_dics_list
	
	stderr('creating lengths distributions at mother_alt_father_ref...')
	mother_alt_father_ref = parents_gt.query(	'(' + maternal_sample_name + ' == 2) & \
							(' + paternal_sample_name + ' == 0) & \
							(indel == 0)').index.values
	mother_alt_father_ref_dics_list = create_length_distributions(mother_alt_father_ref)
	for tup in mother_alt_father_ref_dics_list:
		shared_fragments_dic.update(tup[1]) # alt is shared
		fetal_fragments_dic.update(tup[0]) # ref is fetal
	del mother_alt_father_ref_dics_list

	return shared_fragments_dic, fetal_fragments_dic

def calculate_total_fetal_fraction(shared_frags_dict, fetal_frags_dict):
	total_fetal_fraction = (2 * len(fetal_frags_dict)) / (len(shared_frags_dict) + len(fetal_frags_dict))
	stderr('Total fetal fraction: ', total_fetal_fraction)
	return total_fetal_fraction

def get_all_known_fetal_fragments(parents_gt, fetal_frags_dict):
	'''
	creating a pool of known fetal fragments to add as data to the model
	'''

	stderr('creating a list of known fetal fragments...')
	
	stderr('extracting fetal fragments from mother_ref_father_het positions...')
	mother_ref_father_het = parents_gt.query('(' + maternal_sample_name + ' == 0) & \
							(' + paternal_sample_name + ' == 1) & \
							(indel == 0)').index.values
	mother_ref_father_het_dics_list = create_length_distributions(mother_ref_father_het)
	for tup in mother_ref_father_het_dics_list:
		fetal_frags_dict.update(tup[1])
	del mother_ref_father_het_dics_list
	
	stderr('extracting fetal fragments from mother_alt_father_het positions...')
	mother_alt_father_het = parents_gt.query('(' + maternal_sample_name + ' == 2) & \
							(' + paternal_sample_name + ' == 1) & \
							(indel == 0)').index.values
	mother_alt_father_het_dics_list = create_length_distributions(mother_alt_father_het, fb_by_chrom_dir, fetal = 'ref')
	for tup in mother_alt_father_het_dics_list:
		fetal_frags_dict.update(tup[0])
	del mother_alt_father_het_dics_list

	return fetal_frags_dict