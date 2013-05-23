#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import subprocess
import logging
import os

logger = logging.getLogger('enterMPY')
# This is easily changed.
py_script = 'mPyClust.py'

def mpi_cluster(out_dir,seq_dir, num_procs, sim):
	logger.info('Entering second clustering')
	logger.info('Divvied sequence directory given as ' + seq_dir)
	logger.info(str(num_procs) + " processes requested")

	# Create hc_{sim} parent directory
	parent_dir = os.path.join(out_dir,'hc_{}'.format(str(int(sim*100))))
	try:
		os.makedirs(parent_dir)
	except OSError as e:
		logger.warning('Overwrite will occur on parent directory\n\t'+ str(e))

	passed_args = ['mpiexec',
				   '-n', str(num_procs),
				   'python', py_script, parent_dir, seq_dir, str(sim)]
	logging.debug('Subprocess arguments: {}'.format(passed_args))
	subprocess.call(passed_args,shell=False)

	logging.info('Parent directory for second-cluster: {}'.format(parent_dir))
	return parent_dir