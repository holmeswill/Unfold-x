### Set Endianness!
#endian = '>' # big endian
endian = '<' # little
### Set params
prefix = 'fileout.dat' # the part of the binary files' name before '_Kxxx_Sx'
start_bnd = 10        # first band we want to read
end_band  = 20        # Last band
nkp       = 201       # number of kpoints in band plot  
spin_c    = 1         # spin component to parse (1 or 2)

# Import some usefull stuff

import numpy as np
import struct




for k in range(1,nkp+1):

    filename = prefix + '_K{0:03d}_S{1:01d}'.format(k, spin_c)
    
    with open(filename, mode='rb') as file:
        fileContent = file.read()    
        
        for i in range(start_bnd,end_band):
            pkm, eig = struct.unpack(endian+'dd', fileContent[i*16:i*16+16])
            
            # now do whatever you want!
            print(pkm, eig)
