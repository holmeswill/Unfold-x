# coding: utf-8
import numpy as np
import pylab as plt
import sys

prefix = sys.argv[1]
data = np.fromfile('SC/{}.dat'.format(prefix))
bc = np.loadtxt('BC/band.dat.gnu')

nkpt=44
npts=2001

with open('BC/{}.scf.out'.format(prefix),'r') as f:
    for line in f.readlines():
        if '(alat)' in line:
            BCalat = float(line.split()[-2])
        if 'Fermi' in line:
            EFBC = float(line.split()[-2])
            break
with open('SC/{}.scf.out'.format(prefix),'r') as f:
    for line in f.readlines():
        if '(alat)' in line:
            SCalat = float(line.split()[-2])
        if 'Fermi' in line:
            EFSC = float(line.split()[-2])
            break

alatRatio = BCalat/SCalat

A=data.reshape([nkpt,npts,3])
plt.pcolormesh(A[:,0,0]*alatRatio,A[0,:,1]-EFSC,A[:,:,2].T,shading='auto')
plt.scatter(bc[:,0],bc[:,1]-EFBC,s=10,c='tab:orange')


plt.ylim([0 - EFSC, 20 - EFSC  ])
plt.savefig(sys.argv[1]+'.png')
