from collections import OrderedDict
import vcf
import parse_gt

info_dic = OrderedDict([('MATINFO_FORMAT', {	'num': '.',
						'type': 'String',
						'desc': '"Format of maternal sample ganotyping information"',
						'source': '"parental vcf"'}),
			('MAT_INFO', {		'num': '.',
						'type': 'String',
						'desc': '"Maternal sample ganotyping information"',
						'source': '"parental vcf"'}),
			('PATINFO_FORMAT', {	'num': '.',
						'type': 'String',
						'desc': '"Format of paternal sample ganotyping information"',
						'source': '"parental vcf"'}),
			('PAT_INFO', {		'num': '.',
						'type': 'String',
						'desc': '"Paternal sample ganotyping information"',
						'source': '"parental vcf"'}),
			('PARENTS_QUAL', {	'num': '.',
						'type': 'Float',
						'desc': '"Parental samples QUAL score"',
						'source': '"parental vcf"'})])

reserved_formats = ('GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL')


def print_info_or_format_row(info_or_format, field_id, number, field_type, description, source=False):
	line_list = []
	line_list.append('##' + info_or_format +'=<ID=' + field_id)
	line_list.append('Number=' + number) # int, A, R, G, '.'
	line_list.append('Type=' + field_type)
	if source:
		line_list.append('Source="' + source + '"')
	line_list.append('Description="' + description + '">')

	return ','.join(line_list)
	
def printvcf(x, out_path = False, *args, **kargs):
	if out_path:
		with open(out_path, 'w') as f:
			print(x, file = f, *args, **kargs)
	else:
		print(x, *args, **kargs)

def make_header(cfdna_vcf_reader, parents_vcf_reader, fetal_sample_name, info_dic, reserved_formats, out_path = False):
	
	if cfdna_vcf_reader.metadata['reference'] != parents_vcf_reader.metadata['reference']:
		sys.exit('ERROR: the vcf files are not based on the same reference genome!')
	if cfdna_vcf_reader.contigs != parents_vcf_reader.contigs:
		printerr(	'Warning: cfdna and parental vcf files were aligned to same reference file \
				but their contigs are different')

	# print unique header fields
	printvcf(	'##fileFormat=' + cfdna_vcf_reader.metadata['fileformat'],
			'##fileDate=' + strftime('%Y%m%d'),
			'##source=hoobari',
			'##phasing=none', # phasing is not yet supported
			'##reference=' + cfdna_vcf_reader.metadata['reference'],
			'##commandline="' + ' '.join(sys.argv) + '"',
			sep = '\n',
			out_path = False)
	
	# print contigs header fields
	cfdna_contigs_output = []
	for c in cfdna_vcf_reader.contigs.values():
		cfdna_contigs_output.append('##contig=<ID=' + str(c.id) + ',length=' + str(c.length) + '>')
	printvcf('\n'.join(cfdna_contigs_output), out_path = False)

	# print info header fields
	for k in info_dic:
		print_info_or_format_row(	'INFO',
						k,
						info_dic[k]['num'],
						info_dic[k]['type'],
						info_dic[k]['desc'],
						info_dic[k]['source'])

	# print format header fields
	cfdna_infos_names = [i[0] for i in cfdna_vcf_reader.formats.values()]
	for k in reserved_formats:
		if (k == 'GT') or (k == 'GJ'):
			k_source = 'hoobari'
		else:
			k_source = 'cfdna vcf file'

		print_info_or_format_row(	'FORMAT',
						cfdna_vcf_reader.formats[k].id,
						vcf.Writer.counts[vcf_num_parents_vcf_reader.formats[k].num],
						cfdna_vcf_reader.formats[k].type,
						cfdna_vcf_reader.formats[k].desc,
						k_source)

	# print column names
	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]
	printvcf('\t'.join(vcf_columns), out_path = out_path)

def rec_sample_to_string(rec, sample):
	data = rec.genotype(sample).data
	
	gt = str(data.GT)
	dp = str(data.DP)
	ad = ','.join(str(i) for i in data.AD)
	ro = str(data.RO)
	qr = str(data.QR)
	ao = str(data.AO[0])
	qa = str(data.QA[0])
	gl = ','.join(str(i) for i in data.GL)

	format_and_gt_dic = OrderedDict([	('GT', gt),
						('DP', dp),
						('AD', ad),
						('RO', ro),
						('QR', qr),
						('AO', ao),
						('QA', qa),
						('GL', gl)])

	return format_and_gt_dic

def print_var(variant_name, prediction, phred, pos_info_dic, format_and_gt_dic, out_path = False):

	vcf_list = []
	for row in posteriors_df.itertuples():

		row_list = []
		variant_name = row[0]
		
		# columns 1 - 5
		chrom, pos, ref, alt = parse_gt.uid_to_rec(variant_name)
		row_list += [chrom, pos, '.', ref, alt]
		
		# column 6-7
		row_list += [str(phred), '.']

		# column 8
		info_list = [str(k) + '=' + str(info_dic[k]) for k in info_dic]
		row_list += [';'.join(info_list)]

		# columns 9-10
		format_list = list(format_and_gt_dic.keys())
		fetal_gt_list = list(format_and_gt_dic.values())
		row_list += [':'.join(format_list)] + [':'.join(fetal_gt_list)]

		# merge all to one row string
		variant_row = '\t'.join(row_list)
		
		printvcf(variant_row, out_path = out_path)
