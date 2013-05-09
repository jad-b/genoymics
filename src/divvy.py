#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import os

def get_divvy_dir():
	dir = None
	try:
		os.makedirs('output/divvy_seqs')
	except OSError:
		pass

	dir = 'output/divvy_seqs/'
	return dir



def retrieve_seq_by_id(seq_id,seq_file):
	# print "Retrieving on "+seq_id
	# Iterate through, looking for seq_id
	for line in seq_file:
		if line[0] != '>': # Skip non-headers for efficiency
			continue

		#print "\tChecking "+line.rstrip('\n')
		if seq_id in line:
			# print "\tMatch found!"
			break

	if line == '': # Check for EOF
		return None

	target_seq = line			
	for line in seq_file:
		if line[0] == '>': # Quit on a header
			break		
		target_seq += line 

	seq_file.seek(0) # Reset file position
	# print "\tReturning "+target_seq
	return target_seq



def divvy_otus(divvy_dir, otu_table, seq_file, seqs_per_file ):
	"""
	Assumptions:
	- divvy_dir is real
	- otu_table is a read-only file object for the most-recent otu table, with
	tab-delimited fields
	- seq_file is the base file of sequences
	Reads the assigned sequence IDs for each OTU, and copies them from
	the sequence file to a new file. Attempts to only read a certain number
	of sequences per file, but will not split an OTU.
	"""
	seqs_name_format = seq_file.name.split('/')[-1]
	seqs_name_format = seqs_name_format.split('.')
	base_otu_id = None
	curr_otu_id = 0
	seqs_read = 0
	divvied_seqs = ""

	# print "divvy_otus: seqs_per_file = "+str(seqs_per_file)
	for otu_line in otu_table:

		otu_info = otu_line.rstrip('\n').split('\t')
		curr_otu_id = otu_info[0]
		if base_otu_id is None: 
			base_otu_id = curr_otu_id
		# print "Checking otu "+curr_otu_id+"..."
		for seq_id in otu_info[1:]:
			sequence = retrieve_seq_by_id(seq_id,seq_file)
			if sequence is None:
				raise Exception('divvy: Unable to retrieve '+seq_id)
			divvied_seqs += sequence
			seqs_read += 1
			# print "{} of {}".format(seqs_read,seqs_per_file)
		# Check if we've read enough sequences to write to file
		if seqs_read >= seqs_per_file:
			suffix = "_{}_{}.{}".format(base_otu_id,
										curr_otu_id,
										seqs_name_format[-1])
			out_filename = divvy_dir+seqs_name_format[0]+suffix
			print "divvy: writing {}...".format(out_filename)
			with open(out_filename,'w') as out_file:
				out_file.write(divvied_seqs)
			# print "\tWrite successful"
			base_otu_id = None
			divvied_seqs = ""
			seqs_read = 0
	else: # Write last file
		suffix = "_{}_{}.{}".format(base_otu_id,
										curr_otu_id,
										seqs_name_format[-1])
		out_filename = divvy_dir+seqs_name_format[0]+suffix
		#print "divvy: writing {}...".format(out_filename)
		with open(out_filename,'w') as out_file:
			out_file.write(divvied_seqs)
		# print "\tWrite successful"
		base_otu_id = None
		divvied_seqs = ""
		seqs_read = 0


def sort_seqs_by_otu(otu_table,seq_file,seqs_per_file):
	"""
	Scan the otu table, writing associated sequences to a new
	.fna file. Attempt to assign roughly equivalent numbers of 
	sequences per file, by taking the total number of sequences
	divided by the number of processes to be run.

	Files will be output to the output/divvy_seqs folder.
	Filenames use the sequence filename before the '.fna',
	with the range of OTU numbers suffixed on.
	For example, a divvied file with sequences from OTUs 1 to 122
	and a base sequence filename of seqs_R1_001.fna reads as follow:
	seqs_R1_001_1_122.fna
	"""	

	# Build directory to store these sequences
	divvy_dir = get_divvy_dir()

	# divvy OTUs into separate files
	divvy_otus(divvy_dir, otu_table,seq_file, seqs_per_file)

	#return directory of sorted sequences
	return divvy_dir