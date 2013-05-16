#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import os
import logging

def validate_directory(dir):
	"""
	Check if 'dir' is a valid directory.
	If it exists, return it.
	If it doesn't exist, attempt to make it.
	If an error occurs, raise an exception.
	"""
	try:		
		if os.path.exists(dir):
			return dir
		else:
			os.makedirs(dir)
			return dir
	except OSError as e:
		logging.error('Problem with creation of',dir)
		raise 

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

def write_to_file(out_filename,div_seqs):
	logging.info("divvy: writing {}...".format(out_filename))
	with open(out_filename,'w') as out_file:
		out_file.write(div_seqs)

	return None,"",0

def divvy_otus(divvy_dir, otu_map, seq_file, seqs_per_file ):
	"""
	- divvy_dir is the output directory for divvied sequences
	- otu_dir is the directory containing the otu_map
	- seq_file is the base file of sequences
	- seqs_per_file is the target number of sequences for each file
	Reads the assigned sequence IDs for each OTU, and copies them from
	the sequence file to a new file. Attempts to only read a certain number
	of sequences per file, but will not split an OTU.
	"""
	seqs_dir,seqs_filename = os.path.split(seq_file)
	seqs_name,seq_ext = os.path.splitext(seqs_filename)
	base_otu_id = None
	curr_otu_id = 0
	seqs_read = 0
	file_count = 0
	divvied_seqs = ""

	with open(otu_map,'r') as otus,open(seq_file,'r') as seqs:

		logging.debug("divvy_otus: seqs_per_file = "+str(seqs_per_file))
		for otu_line in otus:

			otu_info = otu_line.rstrip('\n').split('\t')
			curr_otu_id = otu_info[0]
			if base_otu_id is None: 
				base_otu_id = curr_otu_id

			#logging.debug("Checking otu "+curr_otu_id+"...")
			for seq_id in otu_info[1:]:
				sequence = retrieve_seq_by_id(seq_id,seqs)
				if sequence is None:
					raise Exception('divvy: Unable to retrieve '+seq_id)
				divvied_seqs += sequence
				seqs_read += 1
				#logging.debug("{} of {}".format(seqs_read,seqs_per_file))
			# Check if we've read enough sequences to write to file
			if seqs_read >= seqs_per_file:
				logging.debug("Writing at {} over {} sequences".
					format(seqs_read,seqs_per_file))
				out_filename = "{}_{}_{}{}".format(seqs_name, base_otu_id,
					curr_otu_id, seq_ext)
				out_filename = os.path.join(divvy_dir,out_filename)
				# Write and reset all three values
				base_otu_id,divvied_seqs,seqs_read = write_to_file(out_filename,divvied_seqs)
				file_count +=1 

		else: # write last file
			logging.debug("Writing at {} over {} sequences".
					format(seqs_read,seqs_per_file))
			out_filename = "{}_{}_{}{}".format(seqs_name, base_otu_id,
					curr_otu_id, seq_ext)
			out_filename = os.path.join(divvy_dir,out_filename)
			# Write and reset all three values
			base_otu_id,divvied_seqs,seqs_read = write_to_file(out_filename,divvied_seqs)
			file_count +=1 

			return file_count


def sort_seqs_by_otu(otu_map,seq_file,seqs_per_file,
	divvy_output_dir='divvy'):
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
	logging.info("divvying initiated")
	logging.debug("OTU table located in"+otu_map)
	logging.debug("Sequence file located at"+seq_file)
	logging.debug("Target number of sequences per file: "+str(seqs_per_file))
	logging.debug("Target directory for divvied sequences:"+divvy_output_dir)

	# Check directory to store these sequences
	divvy_dir = validate_directory(divvy_output_dir)
	logging.debug('directory validated:'+divvy_dir)

	# divvy OTUs into separate files
	count = divvy_otus(divvy_dir,otu_map,seq_file,seqs_per_file)

	#return directory of sorted sequences
	logging.info("divvying complete")
	return divvy_dir,count
