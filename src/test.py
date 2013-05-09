#!/usr/bin/env python

import subprocess
"""
filename = 'test/tiny454sample.fna'

ret = subprocess.Popen(['grep',">",'-c',filename],
						stdout=subprocess.PIPE,
						shell=False)
print ret.stdout.read().rstrip('\n')
"""
header = ">GKAQVJI01EWIBH rank=0002230 x=1893.5 y=1579.0 length=342"
seq_id = "GKAQVJI01EWIBH"
if seq_id in header:
	print "At least something in this world makes sense."