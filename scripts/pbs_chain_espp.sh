#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=24:00:00
#PBS -N chain_espp_${ll}

module load intel/2014a Boost/1.55.0-intel-2014a-Python-2.7.6 h5py/2.2.1-intel-2014a-Python-2.7.6

cd $PBS_O_WORKDIR
source $HOME/espressopp/ESPRC

make chain_espp PY="mpirun -n 20 python" CHAIN_N=08000 SITES=0001 RUN=${ll} TH=010 CHAIN_STEPS_ESPP=20000 RATE=${RATE} SITES=${SITES}

