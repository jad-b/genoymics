#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import subprocess

# This is easily changed.
py_script = 'mPyClust.py'

def mpi_cluster(seq_dir, num_procs):
	passed_args = ['mpiexec',
				   '-n', str(num_procs),
				   'python',py_script, seq_dir]

	s = subprocess.Popen(passed_args,
						stdout=subprocess.PIPE,
						shell=False)

	return s.stdout.readline()