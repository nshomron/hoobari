from sys import stderr

def printerr(output_to_print, *args, **kargs):
	print(output_to_print, file = stderr, sep = '\t', *args, **kargs)