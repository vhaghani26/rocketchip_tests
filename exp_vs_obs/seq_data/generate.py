import sys
import os

# Path to datasynth program
datasynth = sys.argv[1]

fh = open('to_generate.csv')
param_names = fh.readline().strip().split(',')

for line in fh.readlines():
	endedness, peaktype, testset, name, stacks, rps, padding, stddev, width, length, paired, flank = line.strip().split(',')

	layout = '_'.join([endedness, peaktype, testset])
	
	cmd = f'{datasynth} {name} {stacks} {rps} -o {layout} --padding {padding} --stddev {stddev} --width {width} --length {length} --flank {flank}'
	if paired != 'NA': cmd += f' --paired {paired}'
	try:
		os.makedirs(layout)
		cmd += ' --control input'
		
		os.system(cmd)
		sys.stderr.write(cmd + '\n')	
		os.system(f'mv {layout}/{name}.fa {layout}/genome.fa')
	except:
		cmd += f' --genome {layout}/genome.fa'
		os.system(cmd)
		sys.stderr.write(cmd + '\n')
	os.system(f'gzip {layout}/*.fastq')
