#!/bin/bash

mpirun -np 1 $QE_ROOT/bin/pw.x -inp graphene.scf > graphene.scf.out
mpirun -np 1 $QE_ROOT/bin/pw.x -inp graphene.bands > graphene.bands.out
$QE_ROOT/bin/bands.x < bands.in > bands.out
