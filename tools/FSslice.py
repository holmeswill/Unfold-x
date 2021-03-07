### Set Endianness!
#endian = '>' # big endian
endian = '<' # little
### Set params
unfold_out = 'O_64_imp0.unfold.out2'
prefix     = 'imp-1o4.dat' # the part of the binary files' name before '_Kxxx_Sx'
start_bnd  = 280        # first band we want to read
end_band   = 320        # Last band
nkp        = 400       # number of kpoints in band plot  
spin_c     = 1         # spin component to parse (1 or 2)
E_F        = 12.04     # Fermi energy
### Parameters for plotting the Fermi surface
w          = 0.005      # with of the gaussian for energies close to E_F
kdirection = 2         # draw normal plane w.r.t. 0 -> x, 1 -> y , 2 -> z
kvalue     = 0.0       # value for cut 

# Import some use full stuff

import numpy as np
import struct
from scipy.interpolate import griddata


## Prase unfold.x output

# Parse BZ
# TODO


# Parse K list

k_point_list=[]

with open(unfold_out, mode='r') as f:
    lines = f.readlines()
    for line in lines:
        if 'Doing' in line:
            data = line.split()
            p = np.array( [ float(data[3]) , float(data[4]) , float(data[5][:-1]) ] )# Create a numpy array with point
            k_point_list.append(p)


if len(k_point_list) < nkp:
    print("Error: some k points have no label!")
    exit


positions = []
values    = []

points_coords = [0,1,2]
del points_coords[kdirection]

for k in range(1,nkp+1):
    
    # select slice
    if np.abs(k_point_list[k-1][kdirection]-kvalue) > 0.001:
        continue

    positions.append(k_point_list[k-1][points_coords])
    filename = prefix + '_K{0:03d}_S{1:01d}'.format(k, spin_c)
    
    d = 0.0
    with open(filename, mode='rb') as file:
        fileContent = file.read()    
        
        for i in range(start_bnd,end_band):
            pkm, eig = struct.unpack(endian+'dd', fileContent[i*16:i*16+16])
            
            # now do whatever you want!
            #print(pkm, eig)
            d += pkm * np.exp(-((eig-E_F)**2) / w**2)
            
    values.append(d)
            


grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
FSslice = griddata(np.array(positions), np.array(values), (grid_x, grid_y), method='linear')

import pylab as plt
plt.imshow(FSslice.T, extent=(0,1,0,1), origin='lower')
plt.show()
