Polymerization in Molecular Dynamics for LAMMPS and ESPResSo++
==============================================================

Copyright © 2014-2015 Pierre de Buyl

cg\_md\_polymerization is a repository containing reproducible
computations for Molecular Dynamics simulations of coarse-grained
particles that undergo a polymerization process. Chain growth and step
growth are demonstrated.

This code is written by Pierre de Buyl and is released under the
modified BSD license that can be found in the file LICENSE.

The appropriate citation when using this algorithm is
P. de Buyl and E. Nies, J. Chem. Phys. **142**, 134102 (2015)
[doi:10.1063/1.4916313](http://dx.doi.org/10.1063/1.4916313)
[arXiv:1409.7498](http://arxiv.org/abs/1409.7498). Version 1.0 of this code was used in the publication.

Requirements and usage
----------------------

- [lammps](http://lammps.sandia.gov) The additional fix `bond/create/random` is
  necessary. It is found at
  <https://github.com/pdebuyl/lammps/tree/fbc_random/src/MC>.
- [ESPResSo++](http://www.espresso-pp.de/) version 1.9.
- [Make](http://www.gnu.org/software/make/)
- [HDF5](http://www.hdfgroup.org/HDF5/)
- [Python](https://www.python.org/) with [NumPy](http://www.numpy.org/),
  [matplotlib](http://matplotlib.org/) and [h5py](http://www.h5py.org/)

To reproduce the computations, invoke the make command.

    source /path/to/espressopp/ESPRC
    make chain_lammps LMP="/path/to/lmp"
    make chain_espp
    make epoxy_lammps LMP="/path/to/lmp"
    make epoxy_espp

The seeds that are needed in the simulations, for the random number generators,
are generated from the `/dev/urandom` device of your computer.

Parameters can be given to the simulations:
- `RATE` is the intrisic reaction rate k.
- `TH` is the interval between executions of the polymerization algorithm.
- Depending on the type of simulation, further parameters are available (number
  of particles, functionality of crosslinkers, number of time steps, etc).

Content of the repository
-------------------------

The Python programs for ESPResSo++ are `code/chain_run.py` and
`code/epoxy_run.py` and the input scripts for LAMMPS are `code/in.chain` and
`code/in.epoxy`.

Basic analysis programs are given for the two types of simulations:
`code/analyse_chain.py` and `code/analyse_epoxy.py`.
