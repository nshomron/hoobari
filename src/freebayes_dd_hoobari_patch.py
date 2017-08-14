'''
receives err output pf freebayes -dd as stdin, and for every assessed position, prints out the reads supporting each genotype
'''

from sys import stdin, argv, stderr, exit
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
parser.add_argument("-b", "--bam_file")
parser.add_argument("-t", "--tmp_dir", default = 'tmp_hb')
parser.add_argument("-r", "--region", default = 'region')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'bgzipped vcf of parents, indexed by tabix')
parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--drop_db", action = 'store_true', help = 'override variants database')
args = parser.parse_args()
# ------------------------------

# --------- functions ---------
def get_parental_genotypes(parents_reader, maternal_sample_name, paternal_sample_name, chrom, position):
	n_rec = 0
	for rec in parents_reader.fetch(chrom, int(position) - 1, int(position)):
		maternal_gt = rec.genotype(maternal_sample_name).data.GT
		paternal_gt = rec.genotype(paternal_sample_name).data.GT
		n_rec += 1
		if n_rec > 1:
			exit('more than one parental variant in the same position')
		
	return maternal_gt, paternal_gt

def get_reads_tlen(bam_reader, chrom, position):
	start_pos = max(0, int(position) - 1000)
	end_pos = min(bam_reader.lengths[bam_reader.references.index(chrom)], int(position) + 1000)
	bam_records_at_position = bam_reader.fetch(chrom, start_pos, end_pos) # include a flanking region, since there's local realignment
	tlen_at_position_dic = {}
	for rec in bam_records_at_position:
		tlen_at_position_dic[rec.query_name] = rec.template_length
	return tlen_at_position_dic

def is_fetal_ref(maternal_gt, paternal_gt):
	if maternal_gt == '0/0' and paternal_gt in ('0/1','1/1'):
		return False
	elif maternal_gt == '1/1' and paternal_gt in ('0/0','0/1'):
		return True

def is_fetal_fragment(genotype, ref, alt, fetal_ref = False):
	
	if ((genotype == ref) and fetal_ref) or ((genotype == alt) and not fetal_ref):
		return 1
	else:
		return 0


def get_var_type(vcf_line):
	var_type = re.findall('TYPE=(.*);', vcf_line)[0]
	if var_type == 'snp':
		return 1
	elif var_type == 'mnp':
		return 2
	elif var_type == 'ins':
		return 3
	elif var_type == 'del':
		return 4
	elif var_type == 'complex':
		return 5
	else:
		return 0

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
vardb = db.Variants(args.drop_db, dbpath = os.path.join(args.tmp_dir, 'hoobari.' + str(args.region) + '.db'))

for line in stdin:
	
	if line.startswith('position: '):
		initiate_var = True
		line = line.split()
		var = line[1]

	elif line.startswith('haplo_obs'):
		if initiate_var:
			chrom, position = var.split(':')
			maternal_gt, paternal_gt = get_parental_genotypes(	parents_reader,
										args.maternal_sample_name,
										args.paternal_sample_name,
										chrom,
										position)
			template_lengths_at_position_dic = get_reads_tlen(bam_reader, chrom, position)
			position_list = []
			initiate_var = False
		line = line.rstrip().split('\t')
		geno = line[3]
		read = line[4].split(':')
		qname = ':'.join(read[1:8])
		isize = int(template_lengths_at_position_dic[qname])
		
		position_list.append([geno, math.fabs(isize), qname])

	elif line.startswith('finished position'):
		if not initiate_var:
			var_type = get_var_type(former_line)
			ref, alt = former_line.split('\t')[3:5]

			for l in position_list:
				genotype = l[0]
				is_fetal = is_fetal_fragment(genotype, ref, alt, fetal_ref = is_fetal_ref(maternal_gt, paternal_gt))
				for_ff = use_for_fetal_fraction_calculation(maternal_gt, paternal_gt, var_type, is_fetal)
				l += [is_fetal, var_type, for_ff]

			vardb.insertVariant(chrom.replace('chr',''), int(position), position_list)
			initiate_var = True

	former_line = line


bam_reader.close()
vardb.con.commit()
vardb.con.close()


# TODO: when the db is complete - each qname where is_fetal=1, all appearances of this qname in the db will change to 1
# TODO: problem! needs to know also if fetal only by genotype (not by other fragments) - 
# for_ff_fetal, for_ff_shared


