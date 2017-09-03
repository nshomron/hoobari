from collections import OrderedDict
import vcf
import parse_gt
import sys
from time import strftime
from stderr import *

reserved_formats = ('GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL', 'PG', 'PP')

def print_info_or_format_row(info_or_format, field_id, number, field_type, description, source=False, output_path = False):
	line_list = []
	line_list.append('##' + info_or_format +'=<ID=' + field_id)
	line_list.append('Number=' + str(number)) # int, A, R, G, '.'
	line_list.append('Type=' + field_type)
	line_list.append('Description="' + description + '">')
	if source:
		line_list.append('Source="' + source + '"')

	printvcf(','.join(line_list), out_path = output_path)
	
def printvcf(x, *args, out_path = False, **kargs):
	if out_path:
		with open(out_path, 'a') as f:
			print(x, file = f, *args, **kargs)
	else:
		print(x, *args, **kargs)

def make_header(cfdna_vcf_reader, parents_vcf_reader, input_command, fetal_sample_name, reserved_formats, output_path = False):
	
	if cfdna_vcf_reader.metadata['reference'] != parents_vcf_reader.metadata['reference']:
		printerr('Warning! are the vcf files based on the same reference genome?')
	if cfdna_vcf_reader.contigs != parents_vcf_reader.contigs:
		printerr('Warning! cfdna and parental vcf files have different contigs')

	# print unique header fields
	printvcf(	'##fileformat=' + cfdna_vcf_reader.metadata['fileformat'],
			'##fileDate=' + strftime('%Y%m%d'),
			'##source=hoobari',
			'##phasing=none', # phasing is not yet supported
			'##reference=' + cfdna_vcf_reader.metadata['reference'],
			'##commandline="' + input_command,
			sep = '\n',
			out_path = output_path)
	
	# print contigs header fields
	cfdna_contigs_output = []
	for c in cfdna_vcf_reader.contigs.values():
		cfdna_contigs_output.append('##contig=<ID=' + str(c.id) + ',length=' + str(c.length) + '>')
	printvcf('\n'.join(cfdna_contigs_output), out_path = output_path)

	# TODO: print filter header fields from parents?

	# print format and info header fields	
	printvcf(
	'##FORMAT=<ID=GT,Number=1,Type=String,Source="hoobari",Description="Genotype",Source="hoobari">',
	'##FORMAT=<ID=DP,Number=1,Type=String,Description="Read Depth">',
	'##FORMAT=<ID=AD,Number=R,Type=String,Description="Number of observation for each allele">',
	'##FORMAT=<ID=RO,Number=1,Type=String,Description="Reference allele observation count">',
	'##FORMAT=<ID=QR,Number=1,Type=String,Description="Sum of quality of the reference observations">',
	'##FORMAT=<ID=AO,Number=A,Type=String,Description="Alternate allele observation count">',
	'##FORMAT=<ID=QA,Number=A,Type=String,Description="Sum of quality of the alternate observations">',
	'##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy",Source="hoobari">',
	'##FORMAT=<ID=PG,Number=G,Type=String,Description="P(Genotype), Per-site genotype prior probabilities",Source="hoobari">',
	'##FORMAT=<ID=PP,Number=G,Type=String,Description="P(Posterior), Per-site genotype posterior probabilities",Source="hoobari">',
	'##INFO=<ID=MGT,Number=1,Type=String,Description="Maternal Genotype">',
	'##INFO=<ID=MGQ,Number=1,Type=Float,Description="Maternal Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">',
	'##INFO=<ID=MGL,Number=G,Type=Float,Description="Maternal Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">',
	'##INFO=<ID=MAD,Number=R,Type=Integer,Description="Maternal Number of observation for each allele">',
	'##INFO=<ID=MDP,Number=1,Type=Integer,Description="Maternal Read Depth">',
	'##INFO=<ID=MRO,Number=1,Type=Integer,Description="Maternal Reference allele observation count">',
	'##INFO=<ID=MQR,Number=1,Type=Integer,Description="Maternal Sum of quality of the reference observations">',
	'##INFO=<ID=MAO,Number=A,Type=Integer,Description="Maternal Alternate allele observation count">',
	'##INFO=<ID=MQA,Number=A,Type=Integer,Description="Maternal Sum of quality of the alternate observations">',
	'##INFO=<ID=PGT,Number=1,Type=String,Description="Paternal Genotype">',
	'##INFO=<ID=PGQ,Number=1,Type=Float,Description="Paternal Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">',
	'##INFO=<ID=PGL,Number=G,Type=Float,Description="Paternal Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">',
	'##INFO=<ID=PAD,Number=R,Type=Integer,Description="Paternal Number of observation for each allele">',
	'##INFO=<ID=PDP,Number=1,Type=Integer,Description="Paternal Read Depth">',
	'##INFO=<ID=PRO,Number=1,Type=Integer,Description="Paternal Reference allele observation count">',
	'##INFO=<ID=PQR,Number=1,Type=Integer,Description="Paternal Sum of quality of the reference observations">',
	'##INFO=<ID=PAO,Number=A,Type=Integer,Description="Paternal Alternate allele observation count">',
	'##INFO=<ID=PQA,Number=A,Type=Integer,Description="Paternal Sum of quality of the alternate observations">',
	'##INFO=<ID=MPQ,Number=.,Type=Float,Description="Maternal and Paternal QUAL score from the parental vcf">',
	sep = '\n',
	out_path = output_path)

	# print column names
	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]
	printvcf('\t'.join(vcf_columns), out_path = output_path)

def rec_sample_to_string(rec, sample):
	data = rec.genotype(sample).data
	#print(data)
	
	if data and data.GT != '.':
		format_and_gt_dic = OrderedDict([])
		format_list = rec.FORMAT.split(':')
		for f in format_list:
			idx = format_list.index(f)
			if f in ('AD', 'GL'):
				value = ','.join(str(i) for i in data[idx])
			elif type(data[idx]) is list:
				value = ','.join([str(i) for i in data[idx]])
			else:
				value = str(data[idx])
			
			format_and_gt_dic[f] = value
	else: 
		format_and_gt_dic = '.'
		
	return format_and_gt_dic

def print_var(rec, phred, pos_info_dic, format_and_gt_dic, out_path = False):

	row_list = []

	# columns 1 - 5
	row_list += [rec.CHROM, str(rec.POS), '.', rec.REF, str(rec.ALT[0])]
	
	# column 6-7
	row_list += [str(phred), '.']

	# column 8
	info_list = [str(k) + '=' + str(pos_info_dic[k]) for k in pos_info_dic]
	row_list += [';'.join(info_list)]

	# columns 9-10
	if format_and_gt_dic == '.':
		row_list += [':'.join(reserved_formats), ':'.join(['.'] * len(reserved_formats))]
	else:
		format_list = list(format_and_gt_dic.keys())
		fetal_gt_list = list(format_and_gt_dic.values())
		row_list += [':'.join(format_list)] + [':'.join(fetal_gt_list)]

	# merge all to one row string
	variant_row = '\t'.join(row_list)
	
	printvcf(variant_row, out_path = out_path)

def unsupported_position(rec, out_path = False):
		
		alt = ','.join([str(i) for i in rec.ALT])

		variant_row = [	rec.CHROM,
				str(rec.POS),
				'.',
				rec.REF,
				alt,
				'0',
				'.',
				'.',
				':'.join(reserved_formats),
				':'.join(['.'] * len(reserved_formats)),]
		
		printvcf('\t'.join(variant_row), out_path = out_path)

