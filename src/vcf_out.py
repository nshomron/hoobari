def print_var(variant_row, *args, **kargs):
	if args.vcf_output is not None:
		print(variant_row, file = open(args.vcf_output, 'w'), *args, **kargs)
	else:
		print(variant_row, *args, **kargs)

def make_vcf_header(fetal_sample_name, vcf_output_path):
	
	vcf_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [fetal_sample_name]

def write_var_to_vcf(vcf_output_path):

	posteriors_df['gt'] = posteriors_df.ix[:,0:3].idxmax(axis=1)
	posteriors_df['phred'] = np.where(posteriors_df['gt'] == 0, (-10) * np.log10(1 - posteriors_df[0]), (-10) * np.log10(posteriors_df[0]))

	vcf_list = []
	for row in posteriors_df.itertuples():

		row_list = []
		variant_name = row[0]
		
		# columns 1 - 5
		var_split = re.split(r':|_|/', variant_name)
		row_list += var_split[0:2] + ['.'] + var_split[2:4]
		
		# column 6-7
		row_list += [str('phred')] + ['.']

		# column 8
		info_list = []
		info_list.append('MATINFO_FORMAT=' + '.')
		info_list.append('MAT_INFO=' + '.')
		info_list.append('PAT_FORMAT=' + '.')
		info_list.append('PAT_INFO=' + '.')
		info_list.append('POSTERIORS=' + (','.join(str(p) for p in list(posteriors))))
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
		
		print_var(variant_row)