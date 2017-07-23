import pickle

def pkl_load(file_path):
	print('found a pkl file, loading it...', sep = '')
	with open(file_path, 'rb') as pkl:
		structure = pickle.load(pkl)
	print('done')
	return structure

def pkl_save(obj, file_path):
	with open(file_path, 'wb') as pkl:
		pickle.dump(obj, pkl)