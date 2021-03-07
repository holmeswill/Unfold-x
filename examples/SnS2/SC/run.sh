#!/bin/bash
#SBATCH -J unfold
#SBATCH -o qe.%j.out
#SBATCH -p debug
#SBATCH -n 96 
#SBATCH -N 4
#SBATCH -t 0-00:30:00

NPROC=$SLURM_NPROCS

QEDIR=PATH/TO/qe-6.6/bin
UNFOLDDIR=PATH/TO/qe-6.6/unfold-x/bin

echo $QEDIR

echo $PW

for p in 2 6 8 12 24;
do
NPOOL=$p
PW="mpirun -np $NPROC $QEDIR/pw.x -npool $NPOOL"
UNFOLD1="mpirun -np $NPROC $UNFOLDDIR/unklist.x -npool $NPOOL"
UNFOLD2="mpirun -np $NPROC $UNFOLDDIR/unfold.x -npool $NPOOL"

$PW < scf.in > scf.out
$UNFOLDDIR/unklist.x < unklist.in > unklist.out

cp band.in band.in_$p
cat unklist.out |  sed -n '/K_POINTS t/,/K_POINTS c/p' | sed -e '1d;$d' >> band.in_$p

$PW < band.in_$p > band.out_$p
$UNFOLD2 < unfold.in > unfold.out
mv scf.out scf.out_$p
mv unklist.out unklist.out_$p
mv unfold.out unfold.out_$p
mv SnS2.dat SnS2.dat_$p

rm -r SnS2.save SnS2.wfc* SnS2.xml
done

exit

