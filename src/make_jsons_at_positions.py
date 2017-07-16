import os
import sys
import json
import math


lengths_file = open(sys.argv[1], 'r')
lengths_dic = {}
for line in lengths_file:
	line_split = line.rstrip().split('\t')
	lengths_dic[line_split[0].strip()] = int(line_split[1].strip())
lengths_file.close()

allele_file = open(sys.argv[2], 'r')
snps_dic = {}
for line in allele_file:
	if line.startswith('chr'):
		snp_position = line.rstrip()
		snps_dic[snp_position] = []
	else:
		line_split = line.rstrip().split('\t')
		qname = str(line_split[0])
		geno = str(line_split[1])
		isize = str(lengths_dic[qname])
		snps_dic[snp_position].append((geno, math.fabs(int(isize)), qname))

c = sys.argv[2]

def load_chromosome_json(jsons_dir, chrom):
	# import lengths and qnames from pickle or json file

	print('loading M' + '12148' + '_plasma_with_qnames ' + chrom + ' json or pkl:' )
	chr_length_and_qnames_json = os.path.join(jsons_dir, 'M' + '12148' + '_plasma_with_qnames.chr' + chrom + '.json')
	print('loading ' + chr_length_and_qnames_json + '... ', end = '')
	with open(chr_length_and_qnames_json, 'r') as json_handle:
		chr_length_and_qnames_data = json.load(json_handle)
	print('done')

	return chr_length_and_qnames_data


tmp_dir = os.path.expanduser(sys.argv[1])
out_jsons_dir = os.path.join(tmp_dir, 'jsons', 'chr' + c + '_snps')
os.makedirs(out_jsons_dir, exist_ok=True)

#chromosomes = list(range(1,23)) + ['X']


snps_dic = load_chromosome_json('/groups/nshomron/tomr/projects/cffdna/simulations/lo/lo_from_tmp/12148/wgs/fb_by_chrom', c)

for snp_position in snps_dic:
	output_json = os.path.join(out_jsons_dir, snp_position + '.json')
	with open(output_json, 'w') as fp:
		json.dump(snps_dic[snp_position], fp)

allele_file.close()	
