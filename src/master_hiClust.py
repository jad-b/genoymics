#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

"""
master_hiClust.py coordinates the actions, and is responsible for calling each step. 
This is to provide easy high-level abstraction, as well as to hopefully make each 
piece modular.
"""

#from qiime import pick_otus # from QIIME
import argparse, subprocess			# Standard Python modules
import divvy, enterMPY, consolidate # Custom modules


#Setup cmd line parser
parser = argparse.ArgumentParser(description='Runs hierarchical clustering on a given\
	sequence file')
parser.add_argument('sequences',
					help='The .fna sequence file to cluster')
parser.add_argument('processors',type=int,
					help='The number of processes to run via MPI')
parser.add_argument('-1','--first',default='0.90',
					help='Similarity to cluster for first pass')
parser.add_argument('-2','--second',default='0.97',
					help='Similarity to cluster for second pass')
parser.add_argument('-d','--debug',action='store_true',default=False,
					help='Print debug messages')
parser.add_argument('--skipUclust',action='store_true',default=False,
					help='Skips clustering. Only useful for testing')
parser.add_argument('-sz','--sterilize',type=int,choices=[0,1,2],default=0,
					help='Deletes intermediate files')
args = parser.parse_args()


# variables
sim_first = args.first
first_num = sim_first.split('.')[-1]
sim_second = args.second
second_num = sim_second.split('.')[-1]
method = 'pick_otus.py'
debug = args.debug

# TODO Initiate profiler

# 1) Cluster at 90% 
# Still using subprocess to call QIIME's pick_otus script
if not args.skipUclust:
	uclust_args = [method,
				'-i','{}'.format(args.sequences),
				'-o','output/{}_uclust_otus'.format(first_num),
				'-s',sim_first,
				'-A']
	if debug:
		print uclust_args
	subprocess.call(uclust_args,shell=False)	
	if debug:
		print "First clustering complete"


# 2) Separate sequences for re-clustering

#Get number of sequences in the file
grep_args = ['grep','-c',">",args.sequences]
proc = subprocess.Popen(grep_args,stdout=subprocess.PIPE,shell=False)
num_seqs = int(proc.stdout.read().rstrip('\n'))

# Prepare names
otus_filename = args.sequences.split('/')[-1]
otus_filename_prefix = otus_filename.split('.')[0]
otu_table_path = 'output/{}_uclust_otus/{}_otus.txt'.format(
					sim_first.split('.')[-1],otus_filename_prefix)
# Open files
otu_table = open(otu_table_path,'r')
seqs = open(args.sequences,'r')
print "procs / num_seqs = "+ str(args.processors) + "/" + str(num_seqs)
# Sort sequences by OTU
seq_dir = divvy.sort_seqs_by_otu(otu_table,
								seqs,
								num_seqs/args.processors)
# Close down resources
otu_table.close()
seqs.close()
if debug:
	print "Sequences divvied"

# 3) Re-cluster in parallel
output = enterMPY.mpi_cluster(seq_dir,args.processors).rstrip('\n')
if debug:
	print "Second clustering complete"

# 4) Consolidate sequences by OTU
hiclust_dir = 'hc_{}_otus'.format(sim_second.split('.')[-1])
consol_dir = consolidate.concat_otus(output,
									args.sequences)
if debug:
	print "OTUs consolidated"


if args.sterilize ==  1:
	subprocess.call('rm','-rf','output/divvy_seqs/*',shell=False)
	subprocess.call('rm','-rf','output/*uclust*',shell=False)
elif args.sterilize == 2:
	subprocess.call('rm','-rf','output',shell=False)


