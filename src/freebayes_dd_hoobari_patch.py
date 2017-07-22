'''
receives err output pf freebayes -dd as stdin, and for every assessed position, prints out the reads supporting each genotype
'''

from sys import stdin, argv
import os
import math
import pysam
import argparse
from json_commands import *

# --------- parse args ---------
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--range")
parser.add_argument("-s", "--sample_id")
parser.add_argument("-b", "--bam_file")
parser.add_argument("-t", "--tmp_dir")
args = parser.parse_args()
# ------------------------------

out_dir = os.path.join(args.tmp_dir, 'jsons')
os.makedirs(out_dir, exist_ok=True)
chromosomes = ['chr' + str(i) for i in list(range(1,23)) + ['X']]
[os.makedirs(os.path.join(out_dir, c), exist_ok=True) for c in chromosomes]

bam_reader = pysam.AlignmentFile(os.path.join(args.bam_file), 'rb')

for line in stdin:
	
	if line.startswith('position: '):
		initiate_json = True
		line = line.split()
		var = line[1]
				
	elif line.startswith('haplo_obs'):
		if initiate_json:
			chrom, position = var.split(':')
			bam_records_at_position = bam_reader.fetch(chrom, int(position) - 1000, int(position) + 1000) # take a large flanking area, since there's realignment
			template_lengths_at_position_dic = {}
			for rec in bam_records_at_position:
				template_lengths_at_position_dic[rec.query_name] = rec.template_length
			
			position_file_path = os.path.join(out_dir, chrom, position + '.json')
			position_list = []
			initiate_json = False

		line = line.rstrip().split('\t')
		geno = line[3]
		read = line[4].split(':')
		qname = ':'.join(read[1:8])
		isize = int(template_lengths_at_position_dic[qname])
		position_list.append([geno, math.fabs(isize), qname])

	elif line.startswith('finished position'):
		json_dump(position_list, position_file_path)
		initiate_json = True