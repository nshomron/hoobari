from sys import stderr
from arguments import args as main_args

def printerr(output_to_print, *args, **kargs):
	print(output_to_print, file = stderr, sep = '\t', *args, **kargs)

def printverbose(output_to_print, *args, **kargs):
	if main_args.verbosity:
		print(output_to_print, file = stderr, *args, **kargs)