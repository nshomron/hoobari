'''
receives err output pf freebayes -dd as stdin, and for every assessed position, prints out the reads supporting each genotype
'''

import sys
import os
import re
import math
import pysam
import vcf
import db
import argparse
import pandas as pd

# --------- parse args ---------
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam_file", help = 'The maternal cfDNA bam file\'s path')
parser.add_argument("-t", "--tmp_dir", default = 'tmp_hb')
parser.add_argument("-r", "--region", default = False)
parser.add_argument("-d", "--debug", action = 'store_true', default = False, help = """By default, Freebayes' stderr is used by this patch,
						which can cause a problem when actually trying to debug.
						This flag causes this information to be printed. Please
						note that the data, if saved, ends up in a very large file.""")
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'bgzipped vcf of parents, indexed by tabix')
parser.add_argument("-m", "--m_bam", help = 'maternal bam file')
parser.add_argument("-p", "--p_bam", help = 'paternal bam file')
parser.add_argument("-db", "--db", default = 'hoobari', help = 'db name, or db prefix if hoobari is run per region')
args = parser.parse_args()
# ------------------------------

var_type_dic = {'.': 0, 'snp': 1, 'mnp': 2, 'ins': 3, 'del': 4, 'complex': 5}


# --------- functions ---------
def get_parental_genotypes(parents_reader, parental_samples, chrom, position):
	maternal_sample_name, paternal_sample_name = parental_samples
	n_rec = 0
	for rec in parents_reader.fetch(chrom, int(position) - 1, int(position)):
		maternal_gt = rec.genotype(maternal_sample_name).data.GT
		paternal_gt = rec.genotype(paternal_sample_name).data.GT
		n_rec += 1
		if n_rec > 1:
			sys.exit('more than one parental variant in the same position')

	return maternal_gt, paternal_gt

def get_reads_tlen(bam_reader, chrom, position):
	start_pos = max(0, int(position) - 1000)
	end_pos = min(bam_reader.lengths[bam_reader.references.index(chrom)], int(position) + 1000)
	bam_records_at_position = bam_reader.fetch(chrom, start_pos, end_pos) # include a flanking region, since there's local realignment
	tlen_at_position_dic = {}
	for rec in bam_records_at_position:
		tlen_at_position_dic[rec.query_name] = math.fabs(int(rec.template_length))
	return tlen_at_position_dic

def get_fetal_allele_type(maternal_gt, paternal_gt):
	if maternal_gt == '0/0' and paternal_gt in ('0/1','1/1'):
		return 'alt'
	elif maternal_gt == '1/1' and paternal_gt in ('0/0','0/1'):
		return 'ref'
	else:
		return False

def is_fetal_fragment(genotype, ref, alt, fetal_allele = False):

	if ((genotype == ref) and fetal_allele == 'ref') or ((genotype == alt) and fetal_allele == 'alt'):
		return 1
	elif ((genotype == alt) and fetal_allele == 'ref') or ((genotype == ref) and fetal_allele == 'alt'):
		return 0
	else:
		return None

def get_parental_samples_names(m_bam, p_bam):
	parents_sample_names = []
	for bam_file in (m_bam, p_bam):
		bam_file_reader = pysam.AlignmentFile(os.path.join(bam_file), 'rb')
		sample_name = bam_file_reader.header['RG'][0]['SM']
		parents_sample_names.append(sample_name)
		bam_file_reader.close()
	return parents_sample_names

def use_for_fetal_fraction_calculation(maternal_gt, paternal_gt, var_type, is_fetal):
	var_in_ff_positions = (maternal_gt == '0/0' and paternal_gt == '1/1') or (maternal_gt == '1/1' and paternal_gt == '0/0')
	var_is_snp = var_type == 1

	if var_in_ff_positions and var_is_snp:
		if is_fetal == 1:
			return 1
		elif is_fetal == 0:
			return 2
		else:
			return 0
	else:
		return 0 # not for ff

'''
This patch uses freebayes' algorithm to create a folders' tree: tmp_folder/jsons/chr[1-22,X,Y].
In each chromosome's folder it creates json files named after the position of the variant they represent.

Usage (in bash, using Python 3.5):
freebayes -dd -f [REFERENCE] [OPTIONS] [cfDNA_BAM_FILE] 2>&1 >[OUTPUT] | python freebayes_dd_hoobari_patch.py -b cfDNA_BAM_FILE -t tmp_folder

Explanation:
1) freebayes has to be run with -dd flag, which print more verbose debugging output (and requires "make DEBUG" for installation)
2) Since freebayes' debug information is written to stderr, 2>&1 redirects it to stdout in order to pipe it
3) the tmp_folder created here is the same one you should later use when you run hoobari

1 - var_type - 5 + null
2 - is_fetal - 0/1
3 - ff - 0, 1, 2 (2 + null)

'''

# Initiate variants database
bam_reader = pysam.AlignmentFile(os.path.join(args.bam_file), 'rb')

parents_reader = vcf.Reader(filename = args.parents_vcf)
if args.region:
	dbpath = os.path.join(args.tmp_dir, args.db + '.' + str(args.region) + '.db')
else:
	dbpath = os.path.join(args.tmp_dir, args.db + '.db')
vardb = db.Variants(dbpath = dbpath)

# get parental sample names
parents_sample_names = get_parental_samples_names(args.m_bam, args.p_bam)

# create sample table in the db
vardb.create_samples_table(parents_sample_names)

# write_variant = False
for line in sys.stdin:
	
	if args.debug:
		print(line, file = sys.stderr, end = '')

	
	if line.startswith('#'):
		print(line, end = '')

	elif line.startswith('position: '):
		initiate_var = True
		line = line.split()
		var = line[1]

	elif line.startswith('haplo_obs'):
		if initiate_var:
			chrom, position = var.split(':')
			maternal_gt, paternal_gt = get_parental_genotypes(	parents_reader,
																parents_sample_names,
																chrom,
																position)
			template_lengths_at_position_dic = get_reads_tlen(bam_reader, chrom, position)
			position_list = []
			initiate_var = False
		line = line.rstrip().split('\t')
		geno = line[3]
		read = line[4].split(':')
		qname = ':'.join(read[1:-10])
		isize = template_lengths_at_position_dic[qname]

		position_list.append([geno, isize, qname])
		# print(position_list)

	elif 'TYPE=' in line:

		line_list = line.rstrip().split('\t')
			
		ref, alt = line_list[3:5]
		one_alt_allele = len(alt.split(',')) == 1

		if one_alt_allele:

			var_type_string = line_list[7].split('TYPE=')[1].split(';')[0]
			var_type = var_type_dic[var_type_string]
			
			format_fields = line_list[8].split(':')

			for l in position_list:
				genotype = l[0]	
				is_fetal = is_fetal_fragment(genotype, ref, alt, fetal_allele = get_fetal_allele_type(maternal_gt, paternal_gt))
				for_ff = use_for_fetal_fraction_calculation(maternal_gt,
															paternal_gt,
															var_type,
															is_fetal)
				l += [is_fetal if is_fetal is not None else 0, var_type, for_ff]

			# print(position_list)
			vardb.insertVariant(chrom.replace('chr',''), int(position), position_list)
			print(line, end = '')

vardb.lengthDists()

bam_reader.close()


