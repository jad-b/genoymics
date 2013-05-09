#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import glob, subprocess, sys
from mpi4py import MPI

# Hardcoded parameters
method = 'pick_otus.py'
sim = '0.97'
sim_num = sim.split('.')[1]
out_dir = 'output/hc_'+sim_num+'_otus/'

def cluster(seq_dir):
	"""
	1) Get rank from Communicator
	2) Get sequence file indexed by rank
	3) Cluster on that sequence file
	"""
	# Get rank from Communicator
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	#Get list of divvied sequence files
	dir_list = glob.glob('{}*.fna'.format(seq_dir))
	#print dir_list

	arg_list = ''
	# Get all the names we'll need:
	method_name = method.split('.')[0]
	# regex: (.+/)+.+(_.+)\.(.+)
	otu_suffix = dir_list[rank].split('/')[-1].split('_',1)[-1].split('.')[0]
	out_path = '{path}{sim}_{method}_{otus}/'.format(path=out_dir,
													  sim=sim_num,
													  method=method.split('.')[0],
													  otus=otu_suffix)
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
	if( rank == 0 ):
		print out_dir


def main():
	seq_dir = sys.argv[1] # Get the seq filepath
	cluster(seq_dir) # Cluster it



if __name__ == "__main__":
    main()