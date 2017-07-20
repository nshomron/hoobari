import os
import sys
from json_commands import *
import math
import pysam
import argparse


# --------- parse args ---------
parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_file")
parser.add_argument("-t", "--tmp_dir")
parser.add_argument("-r", "--supporting_reads_file")

args = parser.parse_args()
# ------------------------------

out_dir = os.path.join(args.tmp_dir, 'jsons')
os.makedirs(out_dir, exist_ok=True)
chromosomes = ['chr' + str(i) for i in list(range(1,23)) + ['X']]
[os.makedirs(os.path.join(out_dir, c), exist_ok=True) for c in chromosomes]

bam_reader = pysam.AlignmentFile(os.path.join(args.bam_file), 'rb')

after_first = False
for line in open(args.supporting_reads_file, 'r'):
	if line.startswith('chr'):
		if after_first:
			json_dump(position_list, position_file_path)
		after_first = True
		
		chrom, position = line.rstrip().split(':')
		bam_records_at_position = bam_reader.fetch(chrom, int(position) - 1000, int(position) + 1000) # take a large flanking area, since there's realignment
		template_lengths_at_position_dic = {}
		for rec in bam_records_at_position:
			template_lengths_at_position_dic[rec.query_name] = rec.template_length
		
		position_file_path = os.path.join(out_dir, chrom, position + '.json')
		position_list = []
	else:
		line_split = line.rstrip().split('\t')
		qname = str(line_split[0])
		geno = str(line_split[1])
		isize = int(template_lengths_at_position_dic[qname])
		position_list.append([geno, math.fabs(isize), qname])
# print last one
json_dump(position_list, position_file_path)
