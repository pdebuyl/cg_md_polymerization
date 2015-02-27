#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=12:00:00
#PBS -N chain_espp

module load intel/2014a Boost/1.55.0-intel-2014a-Python-2.7.6 h5py/2.2.1-intel-2014a-Python-2.7.6

cd $PBS_O_WORKDIR
source $HOME/espressopp/ESPRC

NPROC=$( cat $PBS_NODEFILE  |  wc  -l )
make chain_espp PY="mpirun -n ${NPROC} -f ${PBS_NODEFILE} python" CHAIN_N=10000 RUN=${ll} TH=${TH} CHAIN_STEPS_ESPP=${STEPS} RATE=${RATE} SITES=${SITES}

