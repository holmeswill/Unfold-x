# coding: utf-8
import numpy as np
import pylab as plt

data = np.fromfile('SC/Si-unfold.dat')
bc = np.loadtxt('BC/si.dat.gnu')

# quick n dirty
# plt.imshow(np.reshape(data[2::3], (-1, 601)).T , extent=(0,np.max(bc[:,0]),-5,25), aspect=0.1,origin="lower")

# workable, but hard to see 
# plt.scatter(data[0::3],data[1::3],c=data[2::3],s=10)

# Best
npts=601
nkpt=-1 # automatically determined by reshape
A=data.reshape([nkpt,npts,3])
plt.pcolormesh(A[:,0,0],A[0,:,1],A[:,:,2].T,shading='auto')

plt.scatter(bc[:,0],bc[:,1],s=1,c='tab:orange')

plt.ylim([-5,25])
plt.show()
