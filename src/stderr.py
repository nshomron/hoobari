from sys import stderr

def printerr(output_to_print, *args, **kargs):
	print('[pre-processing]', output_to_print, file = stderr, *args, **kargs)