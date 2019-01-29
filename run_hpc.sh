#!/bin/sh
#PBS -r n
#PBS -N par_diffuse_1
##PBS -l nodes=1:ppn=2
#PBS -l nodes=2:ppn=1
#PBS -j oe
#PBS -o par_diffuse.out
#PBS -l walltime=00:5:00

NPROCS=`wc -l < $PBS_NODEFILE` 
echo "NPROCS = " $NPROCS

module add studio
module add mpi/studio
# export PATH=$PATH:/opt/openmpi/1.6.5/sun/bin/mpirun
# change the cd ... below to reach your executable

cd $PBS_O_WORKDIR 

# count=2^0
/opt/openmpi/1.6.5/sun/bin/mpirun -np $NPROCS ./diffusion.exe