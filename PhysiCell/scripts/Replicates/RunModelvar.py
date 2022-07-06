import os
import subprocess

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

os.chdir('../../')

def model(Sample=-1, Replica=-1):
    # Write input for simulation & execute
    callingModel = ['./COVID19', 'output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)


Samples = len(range(0,10,1))
Replicas = 6
numDataPerRank  = int(Samples/size)
mod = Samples%size 
if ( mod != 0): numDataPerRank = numDataPerRank + 1
data = None
if rank == 0:
    data = np.array(list(range(0,10,1)), dtype='d')
    if ( mod != 0):
      add = -1*np.ones((size-mod),dtype=int)
      data = np.concatenate((data,add),axis=None) 
recvbuf = np.empty(numDataPerRank , dtype='d')
comm.Scatter(data, recvbuf, root=0)

for i in range(recvbuf.shape[0]): 
    for idxRep in range(Replicas):
        print('Rank: ',rank, ', recvbuf received: ',recvbuf[i], ', Replica: ', idxRep)
        if ( recvbuf[i] >= 0 ):
          model(recvbuf[i], idxRep)