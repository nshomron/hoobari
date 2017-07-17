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
from json_commands import *
import argparse
# project's
import parse_gt
from stderr import printerr
import vcfuid
import pprogress
import position

# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_id", help = '')
parser.add_argument("-parents", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-tmp", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-no_pkls", "--no_pkls", action = 'store_true', help = "don't make or use pickle files, usually after loading large objects")
parser.add_argument("-v", "--vcf_output", default = 'stdout', help = 'path for vcf output')

args = parser.parse_args()

printerr(args)

# --------- functions ----------

#temp
err_rate = 0.0003

readers = [vcf.Reader(open(cfdna_vcf_path, 'rb')), vcf.Reader(open(parents_vcf_path, 'rb'))]
cfdna_id = readers[0].samples[0]
mother_id = 'M12148W'
father_id = 'H12148W'
reader = vcf.utils.walk_together(readers)

for tup in reader:
	cfdna_rec, parents_rec = tup
	variant_name = vcfuid.rec_to_uid(cfdna_rec)
	
	# calculate priors
	maternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	paternal_gt = parse_gt.str_to_int(parents_rec.genotype(mother_id).data[0])
	priors = position.calculate_priors(maternal_gt, paternal_gt)

	# calculate likelihoods
	likelihoods = position.calculate_likelihoods(variant, args.tmp_dir,	total_fetal_fraction, fetal_fractions_df, err_rate,	lengths = False, origin = False)







def parse_vcf_to_array(vcf_file_path, pkl_file_path):
	stderr('parsing ' + vcf_file_path + '...')
	if os.path.isfile(pkl_file_path):
		vcf_df = pkl_load(pkl_file_path)
	else:
		if os.path.isfile(vcf_file_path):
			index_length = int(subprocess.getoutput(' '.join(['grep','-v','^#',vcf_file_path,'|','wc','-l'])).split()[0])
			progress_index = 0

			vcf_list = []
			with open(vcf_file_path, 'r') as inp:
				reader = vcf.VCFReader(inp)

				samples = reader.samples
				vcf_list.append(['variant_name'] + samples)

				for record in reader:
					variant_name = [vcf_rec_to_var_uid(record)]
					genotypes = []
					for s in range(len(samples)):
						gt = gt_string_to_int(record.genotype(samples[s]).data[0])
						genotypes.append(gt)
					line = variant_name + genotypes

					vcf_list.append(line)
					progress_index = print_progress(progress_index, index_length)
				
				vcf_array = np.array(vcf_list)
				vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
			pkl_save(vcf_df, pkl_file_path)
		else:
			sys.exit('vcf file not found!')
	return vcf_df

def parse_split_vcf_to_array(chr_vcfs_dir_path, files_regex, pkl_file_path):
	stderr('parsing chromosome vcf files from ' + chr_vcfs_dir_path + '...')
	if os.path.isfile(pkl_file_path):
		vcf_df = pkl_load(pkl_file_path)
	else:
		vcf_files = list(filter(re.compile(files_regex).match, os.listdir(chr_vcfs_dir_path)))

		vcf_list = []
		for vcf_file in vcf_files:
			
			stderr('parsing ' + vcf_file + '...')
			vcf_file_path = os.path.join(chr_vcfs_dir_path, vcf_file)
			index_length = int(subprocess.getoutput(' '.join(['grep','-v','^#',vcf_file_path,'|','wc','-l'])).split()[0])

			progress_index = 0
			with open(vcf_file_path, 'r') as inp:
				reader = vcf.VCFReader(inp)

				samples = reader.samples
				vcf_list.append(['variant_name'] + samples)

				for record in reader:
					variant_name = [vcf_rec_to_var_uid(record)]
					genotypes = []
					for s in range(len(samples)):
						gt = gt_string_to_int(record.genotype(samples[s]).data[0])
						genotypes.append(gt)
					line = variant_name + genotypes

					vcf_list.append(line)
					progress_index = print_progress(progress_index, index_length)
			
			vcf_array = np.array(vcf_list)
			vcf_df = pd.DataFrame(data = vcf_array[1:,1:], index = vcf_array[1:,0], columns = vcf_array[0,1:])
		pkl_save(vcf_df, pkl_file_path)
	return vcf_df

def load_chromosome_json(jsons_dir, chrom):
	# import lengths and qnames from pickle or json file

	stderr('loading M' + sample_id + '_plasma_with_qnames ' + chrom + ' json or pkl:' )

	chr_length_and_qnames_pkl = os.path.join(jsons_dir, '..', 'pkl_files', 'M' + sample_id + '_plasma_with_qnames.' + chrom + '.pkl')
	if os.path.isfile(chr_length_and_qnames_pkl):
		chr_length_and_qnames_data = pkl_load(chr_length_and_qnames_pkl)
	else:
		chr_length_and_qnames_json = os.path.join(jsons_dir, 'M' + sample_id + '_plasma_with_qnames.' + chrom + '.json')
		stderr('loading ' + chr_length_and_qnames_json + '... ', end = '')
		with open(chr_length_and_qnames_json, 'r') as json_handle:
			chr_length_and_qnames_data = json.load(json_handle)
		stderr('done')
		pkl_save(chr_length_and_qnames_data, chr_length_and_qnames_pkl)

	return chr_length_and_qnames_data

def create_length_distributions(variant_list, jsons_dir, fetal = None):
	'''
	variant name - form of chr:position_ref/alt
	fetal - choose either the fetal allele is suppose to be the ref (if mother = 1/1 and father = 0/0)
	or the alt (if mother = 0/0 and father = 1/1)
	'''
	
	shared_frags_dic = {}
	fetal_frags_dic = {}
	
	autosomes = ['chr' + str(i) for i in range(1,23)]
	allosomes = ['chrX']
	chromosomes = autosomes + allosomes

	for c in chromosomes:	
		
		chr_variant_list = [var for var in variant_list if var.startswith(c + ':')]
		genotypes_qnames_lengths_tables = load_chromosome_json(jsons_dir,c)

		progress_index = 0
		index_length = len(chr_variant_list)

		for variant_name in chr_variant_list:
			var_split = re.split(r':|_|/', variant_name)
			
			chrom = var_split[0]
			chrom_pos = var_split[0] + ':' + str(var_split[1])
			ref = var_split[2]
			alt = var_split[3]
			if fetal == 'alt':
				shared_allele, fetal_allele = ref, alt
			elif fetal == 'ref':
				shared_allele, fetal_allele = alt, ref
			
			try:
				pos_data = pd.DataFrame(genotypes_qnames_lengths_tables[chrom_pos])
			
				shared_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == shared_allele].iloc[:,1:3]
				for row in shared_lengths_and_qnames_at_position_df.itertuples():
					shared_frags_dic[row[2]] = row[1]
				
				fetal_lengths_and_qnames_at_position_df = pos_data[pos_data[0] == fetal_allele].iloc[:,1:3]
				for row in fetal_lengths_and_qnames_at_position_df.itertuples():
					fetal_frags_dic[row[2]] = row[1]

				progress_index = print_progress(progress_index, index_length)
			
			except:
				progress_index = print_progress(progress_index, index_length)

		del genotypes_qnames_lengths_tables

	return shared_frags_dic, fetal_frags_dic

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

def retrieve_maf_from_ensembl(variant):
	server = "http://grch37.rest.ensembl.org"

	chrom_pos = re.split(r':|_|/', variant)[0:2]
	chrom = chrom_pos[0]
	pos = chrom_pos[1]
	region = chrom + ':' + pos + '-' + pos

	ext = "/overlap/region/human/" + region + "?feature=variation"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	decoded = r.json()
	rsid = decoded[0]['id']

	ext = "/variation/human/" + rsid  + "?pops=1"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	decoded = r.json()
	maf_tuple = ('maf', decoded['MAF'])

	if maf_tuple[1] is None:
		af_list = []
		for pop in decoded['populations']:
			if pop['allele'] == decoded['ancestral_allele']:
				af_list.append(1 - pop['frequency'])
			elif pop['allele'] == decoded['minor_allele']:
				af_list.append(pop['frequency'])

		maf_tuple = ('emaf', np.mean(af_list))
	
	return maf_tuple

def calculate_prior_probabilities(variant_set, parents_gt_pd_array, maternal_gt, pkl_file_path):

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
	stderr('calculating prior probabilities:')
	if os.path.isfile(pkl_file_path):
		priors_df = pkl_load(pkl_file_path)
	else:
		progress_index = 0
		progress_length = len(variant_set)

		maternal_col_idx = parents_gt_pd_array.columns.values[0]
		paternal_col_idx = parents_gt_pd_array.columns.values[1]
		paternal_genotype_dic = parents_gt_pd_array[maternal_col_idx].to_dict()
		paternal_genotype_dic = parents_gt_pd_array[paternal_col_idx].to_dict()

		fetal_genotypes = [0,1,2]
		priors_list = []
		priors_list.append([''] + fetal_genotypes)
		progress_index = 0
		for variant_name in variant_set:
			maternal_gt = maternal_genotype_dic[variant_name]
			paternal_gt = paternal_genotype_dic[variant_name]
			if (maternal_gt in (0,1,2)) and (paternal_gt in (0,1,2)):
				p_maternal_alt = maternal_gt / 2
				p_paternal_alt = paternal_gt / 2
				priors_source = 'parents_vcf'
			else: # get maf or af or ldaf
				maf = retrieve_maf_from_ensembl(variant_name)
				if maf[1] is not None:
					p_maternal_alt = p_paternal_alt = maf
					priors_source = maf[0]
				else:
					priors = 'unknown'
					priors_source = '.'

			priors_at_position = [	(1-p_maternal_alt)*(1-p_paternal_alt),
									p_maternal_alt*(1-p_paternal_alt) + (1-p_maternal_alt)*p_paternal_alt,
									p_maternal_alt*p_paternal_alt]			
			
			for i in range(len(priors_at_position)):
				if priors_at_position[i] == 0:
					priors_at_position[i] = None
				else:
					priors_at_position[i] = np.log(priors_at_position[i])

			progress_index = print_progress(progress_index, progress_length)
			priors_list.append([variant_name] + priors_at_position)

		priors_array = np.asarray(priors_list)
		priors_df = pd.DataFrame(data = priors_array[1:,1:], index = priors_array[1:,0], columns = priors_array[0,1:])
		pkl_save(priors_df, pkl_file_path)
	return priors_df

def calculate_frag_i_with_certain_genotype_given_fetal_genotypes(fetal_genotypes, frag_genotype, ref, alt, f, err_rate):
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
	variant_set,
	jsons_dir,
	total_fetal_fraction,
	err_rate,
	pkl_file_path,
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
	
	stderr('calculating likelihoods... (lengths = ' + str(lengths) + ' and origin = ' + str(origin) + ')')

	if os.path.isfile(pkl_file_path):
		likelihoods_df = pkl_load(pkl_file_path)
	else:
		
		fetal_fractions_df = kwargs.get('ff_df', None)
		known_fetal_qnames_dic = kwargs.get('fetal_qnames', None)

		fetal_genotypes = [0,1,2]
		likelihoods_list = []
		likelihoods_list.append([''] + fetal_genotypes)
		
		autosomes = ['chr' + str(i) for i in range(1,23)]
		allosomes = ['chrX']
		chromosomes = autosomes + allosomes
		
		for c in chromosomes:	
			genotypes_qnames_lengths_tables = load_chromosome_json(jsons_dir,c)
			
			chr_variant_list = [var for var in variant_set if var.startswith(c + ':')]

			progress_index = 0
			progress_length = len(chr_variant_list)
			
			for variant_name in chr_variant_list:
				var_split = re.split(r':|_|/', variant_name)
				chrom = var_split[0]
				chrom_pos = var_split[0] + ':' + str(var_split[1])
				ref = var_split[2]
				alt = var_split[3]

				if genotypes_qnames_lengths_tables[chrom_pos] is not None:
					pos_data = pd.DataFrame(genotypes_qnames_lengths_tables[chrom_pos])

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


						frag_i_likelihood_list = calculate_frag_i_with_certain_genotype_given_fetal_genotypes(fetal_genotypes, frag_genotype, ref, alt, ff, err_rate)
						if frag_i_likelihood_list is not None:
							fragments_likelihoods_list.append(frag_i_likelihood_list)
					
					fragments_likelihoods_df = np.array(fragments_likelihoods_list)
					log_fragments_likelihoods_df = np.log(fragments_likelihoods_df)
					sum_log_fragments_likelihoods_df = log_fragments_likelihoods_df.sum(axis = 0)

					if isinstance(sum_log_fragments_likelihoods_df, np.ndarray):
						variant_line_to_add = [variant_name] + list(sum_log_fragments_likelihoods_df)
						likelihoods_list.append(variant_line_to_add)
						
					progress_index = print_progress(progress_index, progress_length)

			del genotypes_qnames_lengths_tables

		likelihoods_array = np.array(likelihoods_list)
		likelihoods_df = pd.DataFrame(data = likelihoods_array[1:,1:], index = likelihoods_array[1:,0], columns = likelihoods_array[0,1:])
		pkl_save(likelihoods_df, pkl_file_path)

	return likelihoods_df

def calculate_posterior_probabilities(variant_list, priors_df, likelihoods_df, pkl_file_path):
	if os.path.isfile(pkl_file_path):
		all_posteriors_df = pkl_load(pkl_file_path)
	else:
		progress_index = 0
		progress_length = len(variant_list)

		all_posteriors_array = np.empty((len(variant_list),3), dtype = np.float128)
		index_names = []
		i = 0
		for variant_name in variant_list:
			var_priors = np.asarray(priors_df.loc[variant_name], dtype = np.float128)
			var_likelihoods = np.asarray(likelihoods_df.loc[variant_name], dtype = np.float128)
			var_priors_likelihoods = var_priors + var_likelihoods
			var_priors_likelihoods_c = var_priors_likelihoods - np.min(var_priors_likelihoods[~np.isnan(var_priors_likelihoods)])
			exp_var_priors_likelihoods = np.exp(var_priors_likelihoods_c)
			exp_var_priors_likelihoods[np.isnan(exp_var_priors_likelihoods)] = 0
			sum_exp_var_priors_likelihoods = exp_var_priors_likelihoods.sum()
			posteriors_array = np.asarray((exp_var_priors_likelihoods / sum_exp_var_priors_likelihoods), dtype = np.float128)
			index_names.append(variant_name)
			all_posteriors_array[i] = posteriors_array
			i += 1
			progress_index = print_progress(progress_index, progress_length)

		all_posteriors_df = pd.DataFrame(data = all_posteriors_array, index = index_names, columns = [0,1,2], dtype = np.float128)

		pkl_save(all_posteriors_df, pkl_file_path)

	return all_posteriors_df

def gt_int_to_string(gt):
	if gt == 0:
		gt_str = '0/0'
	elif gt == 1:
		gt_str = '0/1'
	elif gt == 2:
		gt_str = '1/1'
	return gt_str

def make_vcf_df(posteriors_df, sample_id, vcf_output_path):

	stderr('Processing the data to create a vcf file...')
	sample_id_hoobari = 'F' + sample_id + 'hb'

	positions_for_vcf = set(posteriors_df.index.values)

	posteriors_df = posteriors_df.loc[positions_for_vcf]
	posteriors_df['gt'] = posteriors_df.ix[:,0:3].idxmax(axis=1)
	posteriors_df['phred'] = np.where(posteriors_df['gt'] == 0, (-10) * np.log10(1 - posteriors_df[0]), (-10) * np.log10(posteriors_df[0]))
	
	index_length = len(posteriors_df)
	progress_index = 0

	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [sample_id_hoobari]
	vcf_list = []
	for row in posteriors_df.itertuples():

		row_list = []
		variant_name = row[0]
		
		# columns 1 - 5
		var_split = re.split(r':|_|/', variant_name)
		row_list += var_split[0:2] + ['.'] + var_split[2:4]
		
		# column 6-7
		row_list += [str(posteriors_df.ix[variant_name, 'phred'])] + ['.']

		# column 8
		info_list = []
		info_list.append('MATINFO_FORMAT=' + '.')
		info_list.append('MAT_INFO=' + '.')
		info_list.append('PAT_FORMAT=' + '.')
		info_list.append('PAT_INFO=' + '.')
		info_list.append('POSTERIORS=' + (','.join(str(p) for p in list(posteriors_df.ix[variant_name,0:3]))))
		row_list += [';'.join(info_list)]

		# columns 9-10
		format_list = []
		fetal_gt_list = []
		
			# GT:DP:AD:RO:QR:AO:QA:GL

		format_list.append('GT')
		fetal_gt_list.append(gt_int_to_string(posteriors_df.ix[variant_name, 'gt']))
		
		row_list += [':'.join(format_list)]
		row_list += [':'.join(fetal_gt_list)]

		vcf_list.append(row_list)
		
		progress_index = print_progress(progress_index, index_length)

	vcf_data_df = pd.DataFrame(data = vcf_list, index = positions_for_vcf, columns = vcf_columns)

	if args.vcf_output == 'stdout':
		vcf_handle = sys.stdout
	else:
		vcf_handle = open(vcf_output_path, 'w')
	stderr('writing output to ' + vcf_output_path + '...')
	vcf_df.to_csv(vcf_handle, sep = '\t', index = False, doublequote=False)
	vcf_handle.close()




# --------- constants ----------
sample_id = args.sample_id
project_dir = os.path.join(os.path.expanduser('~/nshomron_tomr/projects/cffdna/simulations/lo/lo_from_tmp'), sample_id, 'wgs')
parents_vcf_file = os.path.join(project_dir, sample_id + '_parents_wgs.vcf')
fb_by_chrom_dir = os.path.join(project_dir, 'fb_by_chrom')
error_rate_cfdna_vcf = os.path.join(project_dir, sample_id + '_cfdna_error_positions.vcf')

# --------- pickles ----------
if not args.no_pkls:
	pkls = make_pkls(os.path.join(args.tmp_dir, 'pkl_files'))

# --------- main ----------
# parse vcf files
# parents
parents_gt = parse_split_vcf_to_array(project_dir, '.*parents.*chr.*vcf', pkls['parents_gt_pkl'])
parents_gt = parents_gt.drop([var for var in parents_gt.index.values if var.startswith('chrY')], axis=0)
parents_gt = parents_gt.dropna()
# maternal plasma
cfdna_gt = parse_split_vcf_to_array(fb_by_chrom_dir, '.*cfdna.*chr.*vcf', pkls['cfdna_gt_pkl'])
cfdna_gt = cfdna_gt.drop([var for var in cfdna_gt.index.values if var.startswith('chrY')], axis=0)
cfdna_gt = cfdna_gt.dropna()
# intersect to get positions of interest. this shouldn't drop to many positions, since variant calling was done by same coordinates.
positions_to_predict = set(parents_gt.index.values).intersection(cfdna_gt.index.values)
positions_to_predict = set(parents_gt.query('(M' + sample_id + 'W == 1)').index.values).intersection(cfdna_gt.index.values)

# calculations of fragment length distributions, and fetal fraction
mother_ref_father_alt = parents_gt.query('(M' + sample_id + 'W == 0) & (H' + sample_id + 'W == 2)').index.values
mother_alt_father_ref = parents_gt.query('(M' + sample_id + 'W == 2) & (H' + sample_id + 'W == 0)').index.values
stderr(len(mother_ref_father_alt), len(mother_alt_father_ref))

if os.path.isfile(pkls['shared_fragments_dic_pkl']) and os.path.isfile(pkls['fetal_fragments_dic_pkl']):
	stderr('loading shared fragments length distributions...')
	shared_fragments_dic = pkl_load(pkls['shared_fragments_dic_pkl'])
	stderr('done')
	stderr('loading fetal fragments length distributions...')
	fetal_fragments_dic = pkl_load(pkls['fetal_fragments_dic_pkl'])
	stderr('done')
else:
	stderr('creating lengths distributions at mother_ref_father_alt...')
	mother_ref_father_alt_shared_dic, mother_ref_father_alt_fetal_dic = create_length_distributions(mother_ref_father_alt, fb_by_chrom_dir, fetal = 'alt')
	stderr('creating lengths distributions at mother_alt_father_ref...')
	mother_alt_father_ref_shared_dic, mother_alt_father_ref_fetal_dic = create_length_distributions(mother_alt_father_ref, fb_by_chrom_dir, fetal = 'ref')
	shared_fragments_dic = {**mother_ref_father_alt_shared_dic, **mother_alt_father_ref_shared_dic}
	fetal_fragments_dic = {**mother_ref_father_alt_fetal_dic, **mother_alt_father_ref_fetal_dic}
	pkl_save(shared_fragments_dic, pkls['shared_fragments_dic_pkl'])
	pkl_save(fetal_fragments_dic, pkls['fetal_fragments_dic_pkl'])

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
if os.path.isfile(pkls['known_fetal_fragments_dic_pkl']):	
	known_fetal_fragments_dic = pkl_load(pkls['known_fetal_fragments_dic_pkl'])
else:
	stderr('extracting fetal fragments from mother_ref_father_het positions...')
	if os.path.isfile(pkls['mother_ref_father_het_fetal_dic_pkl']):
		mother_ref_father_het_fetal_dic = pkl_load(pkls['mother_ref_father_het_fetal_dic_pkl'])
	else:
		mother_ref_father_het = parents_gt.query('(M' + sample_id + 'W == 0) & (H' + sample_id + 'W == 1)').index.values
		mother_ref_father_het_shared_dic, mother_ref_father_het_fetal_dic = create_length_distributions(mother_ref_father_het, fb_by_chrom_dir, fetal = 'alt')
		del mother_ref_father_het_shared_dic
		pkl_save(mother_ref_father_het_fetal_dic, pkls['mother_ref_father_het_fetal_dic_pkl'])
	
	stderr('extracting fetal fragments from mother_alt_father_het positions...')
	if os.path.isfile(pkls['mother_alt_father_het_fetal_dic_pkl']):
		mother_alt_father_het_fetal_dic = pkl_load(pkls['mother_alt_father_het_fetal_dic_pkl'])
	else:
		mother_alt_father_het = parents_gt.query('(M' + sample_id + 'W == 2) & (H' + sample_id + 'W == 1)').index.values
		mother_alt_father_het_shared_dic, mother_alt_father_het_fetal_dic = create_length_distributions(mother_alt_father_het, fb_by_chrom_dir, fetal = 'ref')
		del mother_alt_father_het_shared_dic
		pkl_save(mother_alt_father_het_fetal_dic, pkls['mother_alt_father_het_fetal_dic_pkl'])
	
	known_fetal_fragments_dic = {**fetal_fragments_dic, **mother_ref_father_het_fetal_dic, **mother_alt_father_het_fetal_dic}
	pkl_save(known_fetal_fragments_dic, pkls['known_fetal_fragments_dic_pkl'])

# make a dictionary that shows the fetal fraction at each fragment length
stderr('making a dictionary that shows the fetal fraction at each fragment length:')
ff_per_length_df = create_fetal_fraction_per_length_df(fetal_fragment_lengths_list, shared_fragment_lengths_list)

# calculate priors
cfdna_priors_df = calculate_prior_probabilities(positions_to_predict, parents_gt, 1, pkls['cfdna_priors_df_pkl'])

# calculate error rate - REMEMBER TO COPY FROM 11189
#error_rate = calculate_error_rate(error_rate_cfdna_vcf, error_rate_pkl)
error_rate = 0.003

# calculate likelihoods and posteriors
# no lengths
cfdna_likelihoods_df = calculate_likelihoods(positions_to_predict, fb_by_chrom_dir, total_fetal_fraction, error_rate, pkls['cfdna_likelihoods_df_pkl'])
predictable_positions = set(cfdna_likelihoods_df.index.values).intersection(cfdna_priors_df.index.values)
stderr('calculating posteriors...')
posterior_probabilities = calculate_posterior_probabilities(predictable_positions, cfdna_priors_df, cfdna_likelihoods_df, pkls['cfdna_posteriors_df_pkl'])

# with lengths, no origin
cfdna_likelihoods_wl_df = calculate_likelihoods(positions_to_predict, fb_by_chrom_dir, total_fetal_fraction, error_rate, pkls['cfdna_likelihoods_wl_df_pkl'],
												lengths = True, ff_df = ff_per_length_df)
stderr('calculating posteriors (with lengths info)...')
posterior_probabilities_wl = calculate_posterior_probabilities(predictable_positions, cfdna_priors_df, cfdna_likelihoods_wl_df, pkls['cfdna_posteriors_wl_df_pkl'])

# with lengths and origin
cfdna_likelihoods_wlo_df = calculate_likelihoods(positions_to_predict, fb_by_chrom_dir, total_fetal_fraction, error_rate, pkls['cfdna_likelihoods_wlo_df_pkl'],
												origin = True, ff_df = ff_per_length_df, fetal_qnames = known_fetal_fragments_dic)
stderr('calculating posteriors (with lengths and origin info)...')
posterior_probabilities_wlo = calculate_posterior_probabilities(predictable_positions, cfdna_priors_df, cfdna_likelihoods_wlo_df, pkls['cfdna_posteriors_wlo_df_pkl'])

# create an output vcf file
make_vcf_df(posterior_probabilities_wlo, sample_id, args.vcf_output)