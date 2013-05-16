#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

"""
master_hiClust.py coordinates the actions, and is responsible for calling each step. 
This is to provide easy high-level abstraction, as well as to hopefully make each 
piece modular.
"""

# Standard Python modules
import argparse 
import subprocess			
import re
import glob
import os
import logging
# Custom modules
import divvy
import enterMPY
import consolidate 


#Setup cmd line parser
parser = argparse.ArgumentParser(
	description=('Runs hierarchical clustering on a given sequence file'))
parser.add_argument('sequences',
					help='The .fna sequence file to cluster')
parser.add_argument('processors',type=int,
					help='The number of processes to run via MPI')
parser.add_argument('-o','--output_directory',help='The top-level output directory',
					default='output',dest='output_directory')
parser.add_argument('-1','--first',default='0.75',type=float,
					help='Similarity to cluster for first pass')
parser.add_argument('-2','--second',default='0.97',type=float,
					help='Similarity to cluster for second pass')
parser.add_argument('-d','--debug',action='store_true',default=False,
					help='Print debug messages')
parser.add_argument('--skipUclust',action='store_true',default=False,
					help='Skips first clustering. Only useful for testing')
parser.add_argument('--keep',action='store_true',default=False,dest='keep',
					help='Deletes intermediate files from this run')
parser.add_argument('-al','-append_log',action='store_true',default=False,
					dest='store_log')
args = parser.parse_args()

if args.debug:
	LEVEL = logging.DEBUG
else:
	LEVEL = logging.INFO
if args.store_log:
	FILEMODE='a' 
else:
	FILEMODE='w'

logging.basicConfig(filename='hc.log',filemode=FILEMODE,
	format='%(levelname)s: %(asctime)s: %(message)s',
	level=LEVEL,
	datefmt='%Y %b %d, %H:%M:%S')
logging.info('master_hiClust.py logging initiated')

# variables
sim_first = args.first
first_num = str(int(sim_first*100))
sim_second = args.second
second_num = str(int(sim_second*100))
num_seqs = 0	# Number of sequences provided
method = 'pick_otus.py'
debug = args.debug


def validate_args():
	if not os.path.exists(args.sequences):
		logging.error("Sequence file does not exist")
		raise Error("Sequence file does not exists")

def retrieve_otu_map_filename(otu_dir):
	otu_map = glob.glob(os.path.join(otu_dir,'*_otus.txt'))
	if otu_map is None:
		logging.error('OTU map filename could not be retrieved')
		raise Exception('Failed to locate otu map from first clustering')

	return otu_map[0]


def count_sequences(seq_file):
	#Get number of sequences in the file
	grep_args = ['grep','-c',">",seq_file]
	proc = subprocess.Popen(grep_args,stdout=subprocess.PIPE,shell=False)
	num = int(proc.stdout.read().rstrip('\n'))
	return num


def setup_dir(_top='output'):
	global num_seqs
	top_dir = _top
	# Build run-specific hiClust directory
	if num_seqs is 0:
		num_seqs = count_sequences(args.sequences)
	hc = 'hc_{}_{}_{}'.format(first_num,second_num,num_seqs)
	hc = os.path.join(top_dir,hc)
	divvy = os.path.join(hc,'divvy')
	logging.debug('output directories:\n'+
		top_dir + '\n' + hc + '\n' + divvy )

	try: 
		if not os.path.exists(divvy):
			# makedirs builds needed intermediate directories
			os.makedirs(divvy)
	except OSError:
		logging.warning('error creating directories')

	return top_dir,hc,divvy


def main():
	validate_args()
	num_seqs = count_sequences(args.sequences)
	top_dir,hc_dir,divvy_dir = setup_dir(args.output_directory)
	uclust_dir = '{}_uclust_otus'.format(first_num)
	uclust_dir = os.path.join(hc_dir,uclust_dir)
	logging.info('uClust directory at '+ uclust_dir)

	# 1) Cluster at 90% 
	if not args.skipUclust:
		uclust_args = [method,
			'-i','{}'.format(args.sequences),
			'-o',uclust_dir,
			'-s',str(sim_first),
			'-A']
		subprocess.call(uclust_args,shell=False)	
	logging.info('First clustering complete')

	# 2) Separate sequences for re-clustering
	# Build OTU map filename from uclust directory
	otu_map = retrieve_otu_map_filename(uclust_dir)
	# Split sequences into files by first-cluster relationships
	(seq_dir, num_files) = divvy.sort_seqs_by_otu(otu_map,args.sequences,
		num_seqs/args.processors,divvy_dir)	

	# 3) Re-cluster in parallel
	# Ultimately we need to run as many processes as we have divvied 
	# sequence files - should be the same, however.
	hc_otu_dir = enterMPY.mpi_cluster(hc_dir,seq_dir,num_files,sim_second)
	logging.info('Second clustering complete')

	# 4) Consolidate sequences by OTU
	final_otu_count = consolidate.concat_otus(hc_otu_dir,hc_dir,args.sequences)
	logging.info('OTUs consolidated')
	logging.info(str (final_otu_count) + " OTUs created")

	print "OTUs",final_otu_count


	if not args.keep:
		rm_args = ['rm','-rf',uclust_dir,divvy_dir,hc_otu_dir]
		logging.info('Deleting intermediate files')
		subprocess.call(rm_args,shell=False)


if __name__ == '__main__':
	main()