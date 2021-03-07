# coding: utf-8
from ase.io import read
from ase.build import make_supercell
import numpy as np
from ase.io.espresso import write_espresso_in
from ase.dft.kpoints import bandpath
from copy import deepcopy
from io import StringIO

# Full path is "GMKGALHA"

BC = read('Zn_mp-79_primitive.cif')
TRMAT=np.array([[2,1,0],
            [-1,1,0],
            [0,0,1]])

SC = make_supercell(BC,TRMAT)
PREFIX='Zn'

# QE Settings
scf_input_data = {
"CONTROL" : {
  'calculation' : 'scf',
  'etot_conv_thr' :   2.0E-05,
  'forc_conv_thr' :   1.0E-04,
  'outdir' : './out/',
  'prefix' : PREFIX,
  'pseudo_dir' : '../pseudo/',
  'verbosity ': 'high'
}, \
'SYSTEM': {
  'degauss' :   1.4699723600E-02,
  'ecutrho' :   3.2E+02,
  'ecutwfc' :   4.E+01,
  'occupations' : 'smearing',
  'smearing' : 'cold'
}, \
'ELECTRONS': {
  'conv_thr' :   4.0E-7,
  'electron_maxstep' : 80,
  'mixing_beta' :   4.0E-01
}
}

pseudopotentials = {'Zn': 'zn_pbe_v1.uspp.F.UPF'}

# Adapt input for band structure calculation
band_input_data = deepcopy(scf_input_data)
band_input_data['CONTROL']['calculation'] = 'bands'

# Generate base cell input files

with open('BC/{}.scf'.format(PREFIX),'w') as f:
    write_espresso_in(f,BC,
                        input_data=scf_input_data,
                        pseudopotentials=pseudopotentials,
                        kspacing=0.05,
                        crystal_coordinates=True)

with open('BC/{}.bands'.format(PREFIX),'w') as f:
    write_espresso_in(f,BC,
                        input_data=band_input_data,
                        pseudopotentials=pseudopotentials,
                        kpts={'path':"GMK",'density':20},
                        crystal_coordinates=True)

with open('BC/{}.pp_band'.format(PREFIX),'w') as f:
    f.write("""
&BANDS
  outdir = './out/'
  prefix = '{}'
  filband='band.dat'
/
""".format(PREFIX)
)

# Generate supercell input files

with open('SC/{}.scf'.format(PREFIX),'w') as f:
    write_espresso_in(f,SC,
                        input_data=scf_input_data,
                        pseudopotentials=pseudopotentials,
                        kspacing=0.05,
                        crystal_coordinates=True)

# Collect reciprocal space path from base cell input
with open('BC/{}.bands'.format(PREFIX),'r') as f:
    ORIG_KPTS = ""
    while True:
        line = f.readline()
        if 'K_POINTS' in line:
            while True:
                l = f.readline()
                if l.strip():
                    ORIG_KPTS += l
                else:
                    break
        if not line:
            break

# Co
with open('SC/{}.scf'.format(PREFIX),'r') as f:
    CELL = ""
    while True:
        line = f.readline()
        if 'CELL_PARAMETERS' in line:
            CELL = line + f.readline() + f.readline() + f.readline()
            break
        if not line:
            break

with open('SC/{}.unklist'.format(PREFIX),'w') as f:
    TRMATTXT = '\n'.join([' '.join([str(y) for y in x]) for x in TRMAT.tolist()])
    f.write(
"""
&input_unklist
    ibrav = 0,
/

TRMAT
{TRMAT}

{CELL_PARAMETERS}

UNKPTS crystal
{KPOINTS}
""".format(TRMAT=TRMATTXT, CELL_PARAMETERS=CELL, KPOINTS=ORIG_KPTS) )

with open('SC/{}.bands'.format(PREFIX),'w') as f:
    write_espresso_in(f,SC,
                        input_data=band_input_data,
                        pseudopotentials=pseudopotentials,
                        kpts=None,
                        crystal_coordinates=True)

with open('SC/{}.bands'.format(PREFIX),'w') as f:
    fio = StringIO()
    write_espresso_in(fio,SC,
                        input_data=band_input_data,
                        pseudopotentials=pseudopotentials,
                        kpts=None,
                        crystal_coordinates=True)
    fio.seek(0)
    for line in fio.readlines():
        if not 'K_POINTS' in line:
            f.write(line)

with open('SC/{}.unfold'.format(PREFIX),'w') as f:
    f.write("""
&inputun
  prefix='{PREFIX}'                  ! prefix of scf calculation
  outdir='./o ut/'                   ! output files directory
  Emin=-20                           ! min energy for spectral function
  Emax=20                            ! max energy for spectral function
  DeltaE=0.02                        ! energy steps from Emin to Emax
  w=0.350                            ! broadening for delta function in spectral function
  filout='./{PREFIX}.dat'            ! output file
  symreduce=.false.                  ! use symmetry to reduce the number of nscf kpoints (probably buggy)
  write_pkm=.false.                  ! output P_Km coefficients and corresponding eigenvalues. WARNING: nspin x nkpt files will be created!
  kpathunit='tpiba'                  ! whether the kpath should be constructed from fractional or Cartesin coordinates.Default: same as UNKPTS units
/

TRMAT
{TRMAT}

UNKPTS crystal
{KPOINTS}
""".format(TRMAT=TRMATTXT, KPOINTS=ORIG_KPTS,PREFIX=PREFIX))


# Final script

with open("run_me.sh", 'w') as f:
    f.write("""#!/bin/bash
if [ -z $QE_ROOT ]; then echo "Please set QE_ROOT"; exit; fi

# Simulation of the base cell
cd BC
$QE_ROOT/bin/pw.x -inp {PREFIX}.scf |tee {PREFIX}.scf.out
$QE_ROOT//bin/pw.x -inp {PREFIX}.bands |tee {PREFIX}.bands.out
$QE_ROOT//bin/bands.x -inp {PREFIX}.pp_band |tee {PREFIX}.pp_band.out

# done with base cell, go back and ...
cd ../SC
# ... do self consistency for supercell
$QE_ROOT//bin/pw.x -inp {PREFIX}.scf |tee {PREFIX}.scf.out

# generate list of k points to be used to generate original band structure
../../../bin/unklist.x < {PREFIX}.unklist |tee {PREFIX}.unklist.out

# generate input for band structure simulation
# (using both crystal and Cartesian coordinates)
cp {PREFIX}.bands {PREFIX}.bands_tpiba
cp {PREFIX}.bands {PREFIX}.bands_crystal
grep -A1000 "K_POINTS crystal" {PREFIX}.unklist.out >> {PREFIX}.bands_crystal
echo "K_POINTS tpiba" >> {PREFIX}.bands_tpiba
awk '/K_POINTS t/{{flag=1; next}} /K_POINTS c/{{flag=0}} flag' {PREFIX}.unklist.out >> {PREFIX}.bands_tpiba

# Now run band structure calculations and unfold
$QE_ROOT//bin/pw.x -inp {PREFIX}.bands_crystal |tee {PREFIX}.bands_crystal.out
../../../bin/unfold.x < {PREFIX}.unfold |tee {PREFIX}.unfold.out

# Plot results
cd .. && python plot.py {PREFIX}

# one may sto here, but we test also the second input
cd SC
$QE_ROOT//bin/pw.x -inp {PREFIX}.bands_tpiba |tee {PREFIX}.bands_tpiba.out
../../../bin/unfold.x < {PREFIX}.unfold |tee {PREFIX}.unfold.out
cd .. && python plot.py {PREFIX}
""".format(PREFIX=PREFIX))
