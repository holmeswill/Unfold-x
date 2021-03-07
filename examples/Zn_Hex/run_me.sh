#!/bin/bash
if [ -z $QE_ROOT ]; then echo "Please set QE_ROOT"; exit; fi

# Simulation of the base cell
cd BC
$QE_ROOT/bin/pw.x -inp Zn.scf |tee Zn.scf.out
$QE_ROOT//bin/pw.x -inp Zn.bands |tee Zn.bands.out
$QE_ROOT//bin/bands.x -inp Zn.pp_band |tee Zn.pp_band.out

# done with base cell, go back and ...
cd ../SC
# ... do self consistency for supercell
$QE_ROOT//bin/pw.x -inp Zn.scf |tee Zn.scf.out

# generate list of k points to be used to generate original band structure
../../../bin/unklist.x < Zn.unklist |tee Zn.unklist.out

# generate input for band structure simulation
# (using both crystal and Cartesian coordinates)
cp Zn.bands Zn.bands_tpiba
cp Zn.bands Zn.bands_crystal
grep -A1000 "K_POINTS crystal" Zn.unklist.out >> Zn.bands_crystal
echo "K_POINTS tpiba" >> Zn.bands_tpiba
awk '/K_POINTS t/{flag=1; next} /K_POINTS c/{flag=0} flag' Zn.unklist.out >> Zn.bands_tpiba

# Now run band structure calculations and unfold
$QE_ROOT//bin/pw.x -inp Zn.bands_crystal |tee Zn.bands_crystal.out
../../../bin/unfold.x < Zn.unfold |tee Zn.unfold.out

# Plot results
cd .. && python plot.py Zn

# one may sto here, but we test also the second input
cd SC
$QE_ROOT//bin/pw.x -inp Zn.bands_tpiba |tee Zn.bands_tpiba.out
../../../bin/unfold.x < Zn.unfold |tee Zn.unfold.out
cd .. && python plot.py Zn
