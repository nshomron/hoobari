from stderr import printerr
import subprocess

def reset():
	progress_index = 0
	return progress_index

def get_file_length(file_path):
	grep = 'grep'
	if file_path.endswith('gz'):
		grep = 'zgrep'
	return int(subprocess.getoutput(' '.join([grep,'-v','^#',file_path,'|','wc','-l'])).split()[0])

def get_length(object):
	return(len(object))

def pprogress(progress_index, length):
	progress_index += 1
	#printerr(str(round(100*((progress_index)/length), 3)) + '%' + '\r', end="")
	
	percents = list(range(0, length + 1, length//100))
	
	if progress_index in percents:
		print(str(progress_index), '/', str(length), '(' + str(int(int(progress_index) / int(length))) + '%)')
	elif progress_index == length:
		print('done')


	return progress_index