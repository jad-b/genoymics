#!usr/bin/env python
"""
mpiPyUClust.py - Python script to impose MPI upon yet another Python script.
Spring 2013
Jeremy Dobbins-Bucklad
Created as part of the Oy* team
"""
from mpi4py import MPI
import os

comm = MPI.COMM_WORLD
pid = None


# initialize MPI
#MPI.Init()
size = comm.Get_size() 
# Get this process's rank (ID) within this process group
rank = comm.Get_rank() 
# Retrieve name of system
name = comm.Get_name()

# Test message
print( "I am {} of {} process on {}\n".format(rank,size,name) )

if( rank == 0 ):
    pid = os.getpid()
    print("PID is {}\n".format(pid))

# Broadcast something as a test
if (rank == 0):
    msg = ['why','am','i','awake']
else:
    msg = None

# Send our message to all processes in our comm group from root
# That's why we don't have to put it in an if-root check
msg = comm.bcast(msg, root=0)

if msg is not None:
    msg_str = " ".join(msg) 
    print "Process {0} says \"{1}\"?\n".format(rank,msg_str)
else:
    print "Process {0} did not get the memo.\n".format(rank)

# Block every process in this comm until they hit this line
#comm.Barrier()

# Terminate the MPI environment
#MPI.Finalize()
 

