import json
from os import path

def json_load(path):
	if path.isfile(path):
		with open(path, 'r') as f:
			json_object = json.load(f)
		return json_object
	else:
		return None