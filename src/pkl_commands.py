import pickle
from stderr import printerr

def pkl_load(file_path):
	printerr('found a pkl file, loading it...')
	with open(file_path, 'rb') as pkl:
		structure = pickle.load(pkl)
	printerr('done')
	return structure

def pkl_save(obj, file_path):
	with open(file_path, 'wb') as pkl:
		pickle.dump(obj, pkl)