#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=6:00:00
#PBS -N epoxy_lammps

module load intel/2014a 

cd $PBS_O_WORKDIR
NPROC=$( cat $PBS_NODEFILE  |  wc  -l )

make chain_growth LMP="mpirun -n ${NPROC} -f ${PBS_NODEFILE} ${HOME}/lammps/src/lmp_ompi_icc" RUN=${ll} FUNC=${FUNC}

