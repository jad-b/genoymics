#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import glob
import subprocess
import sys
import os
from mpi4py import MPI

# Hardcoded parameters
method = 'pick_otus.py'
default_sim = 0.97
sim_num = str(int(default_sim*100))
default_out_dir = 'output/hc_'+sim_num+'_otus/'
debug = False

def cluster(out_dir,seq_dir,sim):
	"""
	1) Get rank from Communicator
	2) Get sequence file indexed by rank
	3) Cluster on that sequence file
	"""
	# Get rank from Communicator
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	#Get list of divvied sequence files
	tmp_dir = os.path.join(seq_dir,'*.fna')
	dir_list = glob.glob(tmp_dir)
	if( rank == 0 and debug ):
		print >> sys.stderr, ('Sequence files to cluster on: '+str(dir_list))

	arg_list = ''
	# Get all the names we'll need:
	if debug:
		print >> sys.stderr, ("rank, filename == {}, {}".format(rank,dir_list[rank]))
	file_path,file_name = os.path.split(dir_list[rank])
	file_basename,file_ext = os.path.splitext(file_name)
	dir_name = 'hc_{sim}_{file}_otus/'.format(sim=sim_num,file=file_basename)
	out_path = os.path.join(out_dir,dir_name)

	if debug: 
		print >> sys.stderr, "Out path is {}".format(out_path)

	#Prepare the arguments for the program
	arg_list = 	['-i','{}'.format(dir_list[rank]),
				'-o',out_path,						
				'-s',sim,
				'-A']

	cmd = [method]
	cmd += arg_list

	# Start uClust instances
	subprocess.call(cmd,shell=False)
	#return out_path
	if( rank == 0 and debug ):
		print out_dir


def main():
	out_dir = sys.argv[1] # Get the output path
	seq_dir = sys.argv[2] # Get the filepath to the divvied sequences
	sim = sys.argv[3] # Get the clustering similarity
	cluster(out_dir,seq_dir,sim) # Cluster it



if __name__ == "__main__":
    main()