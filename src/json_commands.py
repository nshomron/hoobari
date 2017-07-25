import json
from os import path

def json_load(file_path):
	if path.isfile(file_path):
		with open(file_path, 'r') as f:
			json_object = json.load(f)
		return json_object
	else:
		return None

def json_dump(json_object, file_path):
	with open(file_path, 'w') as f:
		json.dump(json_object,f)