
def str_to_int(str_gt):
	if str_gt in ('0/0', '0/1', '1/1'):
		gt = int(str_gt[0]) + int(str_gt[2])
	elif str_gt == '.':
		gt = None
	else:
		gt = 'unsupported'

	return gt

def int_to_str(gt):
	if gt == 0:
		string = '0/0'
	elif gt == 1:
		string = '0/1'
	elif gt == 2:
		string = '1/1'
	else:
		string = '.'
	return string
