#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=8:00:00
#PBS -N chain_lammps

module load intel/2014a 

cd $PBS_O_WORKDIR
NPROC=$( cat $PBS_NODEFILE  |  wc  -l )

make chain_lammps LMP="mpirun -n ${NPROC} -f ${PBS_NODEFILE} ${HOME}/lammps/src/lmp_ompi_icc" CHAIN_N=10000 SITES=${SITES} RUN=${ll} TH=${TH} CHAIN_STEPS_LAMMPS=${MEGA}000000 RATE=${RATE} LAMMPS_CHAIN_FILE=${LCF}

