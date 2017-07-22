import json
from os import path

def json_load(path):
	if path.isfile(path):
		with open(path, 'r') as f:
			json_object = json.load(f)
		return json_object
	else:
		return None

def json_dump(object, path):
	with open(path, 'w') as f:
		json.dump(object,f)