Polymerization in Molecular Dynamics
====================================

Copyright © 2015 Pierre de Buyl

cg\_md\_polymerization is a repository containing reproducible
computations for Molecular Dynamics simulations of coarse-grained
particles that undergo a polymerization process. Chain growth and step
growth are demonstrated.

This code is written by Pierre de Buyl and is released under the
modified BSD license that can be found in the file LICENSE.

Requirements and usage
----------------------

- [lammps](http://lammps.sandia.gov)
- [Make](http://www.gnu.org/software/make/)
- [bash](https://www.gnu.org/software/bash/)
- [HDF5](http://www.hdfgroup.org/HDF5/)
- [Python](https://www.python.org/) with [NumPy](http://www.numpy.org/),
  [matplotlib](http://matplotlib.org/) and [h5py](http://www.h5py.org/)

Apart from lammps that requires a separate installation, the other packages can
be installed under Debian with

    apt-get install make bash libhdf5-dev hdf5-tools python python-numpy python-matplotlib python-h5py

The custom dump style for [H5MD](http://nongnu.org/h5md/) is available at
<https://github.com/pdebuyl/lammps> (in the branch `start_dump_h5md` and is
needed to take advantage of the provided analysis tools.

To reproduce the test computation, invoke the make command.

    make step_growth

This runs a simulation of step growth where each particle can link to a maximum
of 8 other particles.

If a lammps executable named `lmp_mpi` is not found in your PATH environment
variable, you may specify it on the command-line

    make step_growth LMP=/path/to/lammps

Specific parameters to the Lennard-Jones system can be appended to the
command-line as `PROB` (probability of making a bond for each candidate pair).
Seeds for the initial velocities and for the fix `bond/create` are generated
from the `/dev/urandom` device of your computer.

