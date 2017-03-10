# The script here is involved in getting base-level mutation rates in the function adjust_mutation_rate_window_to_base
# [Input] is a bed file and each row is a genomic interval (e.g., 50bp) and with a uniq name
# [mutation_rate_file_path] is where the mutation rate bw file is. (e.g., /media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw)

from __future__ import print_function
import sys
import re
import subprocess


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if sys.argv[1] == '--help':
	print('170126_get_baseline_mutation_rate.py [input] [mutation_rate_file_path]')
	sys.exit()

file_in = sys.argv[1]
tab_file = sys.argv[2]
header = sys.argv[3]
fh = open(file_in, 'r')
fh = fh.readlines()
#fh.pop(0) # don't know why but she removes the first line in the code
counter = 1

header = open(header, 'r')
header = header.readlines()
header = re.sub('\n', '', header[0])
header = header.split('\t')
nas = [ 'NA' for i in range(len(header) - 2) ]
# header = 'index\tdpsi_max_tissue\tdpsi_zscore\tgene\tstrand\ttranscript\texon_number\tlocation\tcds_type\tss_dist\tcommonSNP_rs#'
# nas = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
nas = '\t'.join(nas)
header = '\t'.join(header)
print(header)
for i in fh:
	i = re.sub('\n', '', i) # remove the \n at the last position of each line
	i = i.split('\t') # split i into a vector
	chrm = i[0]
	site = i[1]
	ref = i[3]
	alt = i[4]
	idx = i[2]
	cmd = 'tabix ' + tab_file + ' ' + chrm + ':' + site + '-' + site
	# ret = subprocess.call(cmd, shell = True)
	proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
	(ret, err) = proc.communicate()	
	if str(ret) == 'b\'\'':
		eprint('no return at line' + str(counter) + ' in ' + file_in + '(index = ' + idx + ')')
		tmp = idx + '\t'  + nas + '\tout_of_region'
		print(tmp)
	else:
		ret = ret.splitlines()
		flag = 0
		for j in ret:
			#j = j.split('\t')
			j = str(j).split('\\t')
			ref_spidex = j[2]
			alt_spidex = j[3]
			if ref != ref_spidex:
				eprint('something wrong at line' + str(counter) + ' in ' + file_in + '(index = ' + idx + '): ' + 'ref is ' + ref + ' but ref_spidex is ' + ref_spidex)
				tmp = idx + '\t'  + nas + '\twrong_ref'
				print(tmp)
				flag = 1
				break
			else:
				if alt == alt_spidex:
					re_spidex = j[4:]
					re_spidex.insert(0, idx)
					eprint('hit a mutation at line' + str(counter) + ' in ' + file_in + '(index = ' + idx + ')')
					print('\t'.join(re_spidex))
					flag = 1
					break
		if flag == 0:
			eprint('this mutation is not SNV at line' + str(counter) + ' in ' + file_in + '(index = ' + idx + ')')
			tmp = idx + '\t'  + nas + '\tnot_SNV'
			print(tmp)
	counter += 1
