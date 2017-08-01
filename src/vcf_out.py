from collections import OrderedDict
import vcf
import parse_gt
import sys
from time import strftime
from stderr import *

info_dic = OrderedDict([('PARENTS_FORMAT', {	'num': '.',
						'type': 'String',
						'desc': 'Format of parental sample ganotyping information',
						'source': 'parental vcf'}),
			('MAT_INFO', {		'num': '.',
						'type': 'String',
						'desc': 'Maternal sample ganotyping information',
						'source': 'parental vcf'}),
			('PAT_INFO', {		'num': '.',
						'type': 'String',
						'desc': 'Paternal sample ganotyping information',
						'source': 'parental vcf'}),
			('PARENTS_QUAL', {	'num': '.',
						'type': 'Float',
						'desc': 'Parental samples QUAL score',
						'source': 'parental vcf'}),
			('PROB_SOURCE', {	'num': '.',
						'type': 'String',
						'desc': 'Whether the probabilities for the calculations are only the posteriors, only the likelihoods, or joint',
						'source': 'hoobari'})])

reserved_formats = ('GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GJ')

hoobari_formats_dic = OrderedDict([('GJ', {	'num': 'G',
						'type': 'Float',
						'desc': 'Genotype Joint probabilities, log scaled. This is the Bayesian numerator (likelihood*prior) for each possible genotype (0/0, 0/1, 1/1)',
						'source': 'hoobari'})])

def print_info_or_format_row(info_or_format, field_id, number, field_type, description, source=False, output_path = False):
	line_list = []
	line_list.append('##' + info_or_format +'=<ID=' + field_id)
	line_list.append('Number=' + str(number)) # int, A, R, G, '.'
	line_list.append('Type=' + field_type)
	if source:
		line_list.append('Source="' + source + '"')
	line_list.append('Description="' + description + '">')

	printvcf(','.join(line_list), out_path = output_path)
	
def printvcf(x, *args, out_path = False, **kargs):
	if out_path:
		with open(out_path, 'a') as f:
			print(x, file = f, *args, **kargs)
	else:
		print(x, *args, **kargs)

def make_header(cfdna_vcf_reader, parents_vcf_reader, input_command, fetal_sample_name, info_dic, reserved_formats, output_path = False):
	
	if cfdna_vcf_reader.metadata['reference'] != parents_vcf_reader.metadata['reference']:
		printerr('Warning! are the vcf files based on the same reference genome?')
		printerr('cfDNA:', str(cfdna_vcf_reader.metadata['reference']))
		printerr('parental:', str(parents_vcf_reader.metadata['reference']))
	if cfdna_vcf_reader.contigs != parents_vcf_reader.contigs:
		printerr('Warning! cfdna and parental vcf files have different contigs')
		printerr('cfDNA:', str(cfdna_vcf_reader.contigs))
		printerr('parental:', str(parents_vcf_reader.contigs))

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


	# TODO: print filter header fields from parents

	# print format header fields
	cfdna_infos_names = [i[0] for i in cfdna_vcf_reader.formats.values()]
	for k in reserved_formats:
		
		k_class, k_id = 'FORMAT', k

		if k in cfdna_infos_names:
			
			if cfdna_vcf_reader.formats[k].num in vcf.Writer.counts:
				k_num = vcf.Writer.counts[cfdna_vcf_reader.formats[k].num]
			else:
				k_num = cfdna_vcf_reader.formats[k].num

			k_type = str(cfdna_vcf_reader.formats['GT'].type)
			k_desc = cfdna_vcf_reader.formats[k].desc
			k_source = 'cfdna vcf file'

			if k == 'GT':
				k_source = 'hoobari'
				
		else:
			k_num = hoobari_formats_dic[k]['num']
			k_type = hoobari_formats_dic[k]['type']
			k_desc = hoobari_formats_dic[k]['desc']
			k_source = 'hoobari'
		
		print_info_or_format_row(k_class, k_id, k_num, k_type, k_desc, source = k_source, output_path = output_path)

	# print info header fields
	for k in info_dic:
		print_info_or_format_row(	'INFO',
						k,
						info_dic[k]['num'],
						info_dic[k]['type'],
						info_dic[k]['desc'],
						info_dic[k]['source'],
						output_path = output_path)

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
		row_list += ['.', '.']
	else:
		format_list = list(format_and_gt_dic.keys())
		fetal_gt_list = list(format_and_gt_dic.values())
		row_list += [':'.join(format_list)] + [':'.join(fetal_gt_list)]

	# merge all to one row string
	variant_row = '\t'.join(row_list)
	
	printvcf(variant_row, out_path = out_path)

def unsupported_position(rec, out_path = False):
		variant_row = [	rec.CHROM,
				str(rec.POS),
				'.',
				rec.REF,
				rec.ALT,
				'.',
				'.',
				'MATINFO_FORMAT=.;MAT_INFO=.;PATINFO_FORMAT=.;PAT_INFO=.;PARENTS_QUAL=.',
				'.',
				'.']
		
		printvcf(variant_row, out_path = out_path)

