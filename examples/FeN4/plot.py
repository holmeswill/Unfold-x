# coding: utf-8
import numpy as np
import pylab as plt

data = np.fromfile('SC/FeN.dat')
bc = np.loadtxt('BC/bd.dat.gnu')

# quick n dirty
#plt.imshow(np.reshape(data[2::3], (-1, 3001)).T , extent=(0,np.max(bc[:,0]),0,30), aspect=0.1,origin="lower")

# much better
npts=3001
nkpt=-1 # automatically determined by reshape
A=data.reshape([nkpt,npts,3])
x = A[:,0,0]

# remove jumps in k path
avg_step = np.average(np.diff(A[:,0,0]))
x[1:] -= np.cumsum(np.where(np.diff(A[:,0,0]) > 2*avg_step,np.diff(A[:,0,0]),0))

plt.pcolormesh(x,A[0,:,1],A[:,:,2].T,shading='auto')

plt.scatter(bc[:,0],bc[:,1],s=1,c='tab:orange')


plt.ylim([0,30])
plt.savefig('FeN4.png')
plt.show()
