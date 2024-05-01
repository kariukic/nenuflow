#! /usr/bin/env python3
from argparse import ArgumentParser, FileType
import numpy as np
import sys
import os

mama = "/net/node115/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01.MS"

nodes = np.arange(106, 114)
timesteps = np.arange(0, 6720, 840)

n_tsteps=840


for i, (n, t) in enumerate(zip(nodes, timesteps)):
    node_dir = f"/net/node{n}/data/users/lofareor/chege/nenufar/obs/L2_BP"
    os.system(f"rm -r {node_dir}")
    os.makedirs(node_dir, exist_ok=True)
    msout=f"{node_dir}/SW01_56MinChunk.MS"
    tstart = (0 if i==0 else t+1)
    
    comm = f"DP3 msin={mama} msin.datacolumn=DATA msin.starttimeslot={tstart} msin.ntimes={n_tsteps} msout={msout} msout.datacolumn=DATA steps=[]"

    print(i+1, n, tstart, msout, comm)

    os.system(comm)