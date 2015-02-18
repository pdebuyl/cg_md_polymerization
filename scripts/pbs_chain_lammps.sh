#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=24:00:00
#PBS -N chain_lammps_$${ll}

module load intel/2014a 

cd $PBS_O_WORKDIR

make chain_growth LMP="mpirun -n 20 $HOME/lammps/src/lmp_ompi_icc" CHAIN_N=08000 SITES=0030 RUN=${ll} TH=010 CHAIN_STEPS_LAMMPS=8000000 RATE=${RATE} SITES=${SITES}

