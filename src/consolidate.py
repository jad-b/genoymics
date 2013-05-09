#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import glob

def concat_otus(otu_dir,seq_filename):
	"""
	Compiles the n OTU tables written into one file and writes it 
	to hc_{seqs}_otus.txt
	"""
	seq_name = seq_filename.split('/')[-1].split('.')[0]
	otu_filename = 'hc_{seqs}_otus.txt'.format(seqs=seq_name)
	print "consolidate: " + otu_dir+otu_filename
	with open(otu_dir+otu_filename,'w') as otu_file:
		otu_num = 0
		#print "\t inside 'with'"
		# get all otu filepaths
		dir_list = glob.glob('{}*/*.txt'.format(otu_dir))
		#print "consolidate.concat_otus: " + str(dir_list)
		# Iterate over files
		for f_path in dir_list:
			with open(f_path,'r') as curr_f:
				for line in curr_f:
					otu_seqs = line.split('\t',1)[1]
					to_write = str(otu_num)+' '+ otu_seqs
					otu_file.write(to_write)
					otu_num += 1
	