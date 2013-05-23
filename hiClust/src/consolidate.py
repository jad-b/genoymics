#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import glob
import logging
import os

logger = logging.getLogger('Consolidate')

def concat_otus(otu_dir,output_dir,seq_file):
	"""
	Compiles the n OTU tables into one file and writes it 
	to hc_{seq filename}_otus.txt
	Returns a tuple of the consolidated OTU map and the number 
	of OTUs written.
	"""
	logger.info("Consolidation initiated")

	# Build filepath
	logger.debug("seq_file is " + seq_file)
	seqs_path,seqs_filename = os.path.split(seq_file)
	logger.debug(seqs_path + ", " + seqs_filename)
	seqs_basename,seqs_ext = os.path.splitext(seqs_filename)
	logger.debug(seqs_basename + ", " + seqs_ext )
	otu_filename = 'hc_{seqs}_otus.txt'.format(seqs=seqs_basename)
	otu_dest = os.path.join(output_dir,otu_filename)
	logger.info("Consolidating output at " + otu_dest)

	with open(otu_dest,'w') as otu_file:
		otu_num = 0
		# get all otu filepaths		
		tmp_dir = os.path.join(otu_dir,"*/*.txt")
		dir_list = glob.glob(tmp_dir)
		logger.debug("Consolidating these files:\n" + str(dir_list))
		logger.debug("Consolidating files under " + otu_dir )

		for f_path in dir_list:
			with open(f_path,'r') as curr_f:
				logger.info("Consolidating " + curr_f.name)
				files_otus = 0
				for line in curr_f:
					otu_seqs = line.split('\t',1)[1]
					to_write = str(otu_num)+' '+ otu_seqs
					otu_file.write(to_write)
					otu_num += 1
					files_otus += 1
				logging.debug('Contains {} OTUs'.format(files_otus))
	logger.info("Wrote {} OTUs to {}".format(otu_num,otu_dest))
	return otu_dest,otu_num 

	