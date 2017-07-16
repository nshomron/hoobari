from re import split as regsplit

def rec_to_uid(rec):
	uid = rec.CHROM + ':' + str(rec.POS) + '_' + rec.REF + '/' + str(rec.ALT[0])
	return uid

def uid_to_rec(uid):
	rec = regsplit(r':|_|/', uid)
	return rec