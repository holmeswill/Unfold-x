#!/bin/bash

mpirun -np 1 $QE_ROOT/bin/pw.x -inp Si.scf > Si.scf.out

../../../bin/unklist.x < Si.unklist > Si.unklist.out

cp Si.bnd.template Si.bnd
cat Si.unklist.out |  sed -n '/K_POINTS t/,/K_POINTS c/p' | sed -e '1d;$d' >> Si.bnd

mpirun -np 1 $QE_ROOT/bin/pw.x -inp Si.bnd > Si.bnd.out
../../../bin/unfold.x < Si.unfold > Si.unfold.out
