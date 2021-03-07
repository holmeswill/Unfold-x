# coding: utf-8
import numpy as np
import pylab as plt

data = np.fromfile('SC/graphene.dat')
bc = np.loadtxt('BC/bands.dat.gnu')


plt.imshow(np.reshape(data[2::3], (-1, 3001)).T , extent=(0,np.max(bc[:,0]),-22,8), aspect=0.1,origin="lower")
plt.scatter(bc[:,0],bc[:,1],s=1,c='tab:orange')


plt.ylim([-22,8])
plt.show()
