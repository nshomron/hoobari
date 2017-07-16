from stderr import printerr
import subprocess

def reset():
	progress_index = 0
	return progress_index

def get_file_length(file_path):
	return int(subprocess.getoutput(' '.join(['grep','-v','^#',file_path,'|','wc','-l'])).split()[0])

def get_length(object):
	return(len(object))

def pprogress(progress_index, length):
	progress_index += 1
	printerr(str(round(100*((progress_index)/length), 3)) + '%' + '\r', end="")

	if progress_index == length:
		printerr('done')

	return progress_index