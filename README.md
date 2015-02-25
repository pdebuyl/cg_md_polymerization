Polymerization in Molecular Dynamics for LAMMPS and ESPResSo++
==============================================================

Copyright © 2014-2015 Pierre de Buyl

cg\_md\_polymerization is a repository containing reproducible
computations for Molecular Dynamics simulations of coarse-grained
particles that undergo a polymerization process. Chain growth and step
growth are demonstrated.

This code is written by Pierre de Buyl and is released under the
modified BSD license that can be found in the file LICENSE.

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
