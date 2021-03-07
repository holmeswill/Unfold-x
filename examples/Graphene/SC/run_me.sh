#!/bin/bash

mpirun -np 1 $QE_ROOT/bin/pw.x -inp graphene.scf > graphene.scf.out

../../../bin/unklist.x < graphene.unklist > graphene.unklist.out

cp graphene.bands.template graphene.bands
cat graphene.unklist.out |  sed -n '/K_POINTS t/,/K_POINTS c/p' | sed -e '1d;$d' >> graphene.bands

mpirun -np 1 $QE_ROOT/bin/pw.x -inp graphene.bands > graphene.bands.out
../../../bin/unfold.x < graphene.unfold > graphene.unfold.out
