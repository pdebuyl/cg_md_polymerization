"""Chain growth Molecular Dynamics
===============================

This program performs the following kind of simulation:
* Chain growth of a single polymer chain.
* Chain growth of multiple polymer chains.
* Arrested chain growth for a single chain: using the `--stop-at` command-line
  argument, one can restrict the number of bonds that are generated.

All runs start by a warmup stage that is not recorded to the output file.

The trajectory is dumped during polymerization, except for arrested chain growth
where it is dumped for subsequents time steps.

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int, help='Number of monomers')
parser.add_argument('--dt', type=float, help='Timestep.', default=0.01)
parser.add_argument('--warmup-steps', type=int, help='Number of steps in warmup loop.', default=100)
parser.add_argument('--warmup-loops', type=int, help='Number of iterations in warmup loop.', default=100)
parser.add_argument('--steps', type=int, help='Number of steps in main loop.', default=100)
parser.add_argument('--loops', type=int, help='Number of iterations in main loop.', default=100)
parser.add_argument('--density', type=float, help='Number density.', default=0.8)
parser.add_argument('--sites', type=int, help='Number of initial active sites.', default=1)
parser.add_argument('--rate', type=float, help='Reaction rate.', default=0.1)
parser.add_argument('--interval', type=int, help='Steps between reactions.', default=100)
parser.add_argument('--file', type=str, help='Filename for H5MD output.', default=None)
parser.add_argument('--dump-interval', type=int, help='Interval at which to dump the configuration.', default=20)
parser.add_argument('--stop-at', type=int, help='Length at which a single chain should stop.', default=-1)
parser.add_argument('--seed', type=int, help='Seed for the rng.')
args = parser.parse_args()

if args.stop_at>0:
    assert args.sites==1
    assert args.interval>=args.steps and args.interval%args.steps==0
elif args.stop_at==0:
    assert args.sites==0

import espresso
import chain_setup
import random

system, integrator, LJCapped, verletList, thermostat, num_particles = chain_setup.monomer_system(args.N, density=args.density, seed=args.seed)

for i in random.sample(range(args.N), args.sites):
    system.storage.modifyParticle(i, 'state', 2)

integrator.dt = args.dt

epsilon_start = 0.1
epsilon = 1.0

espresso.tools.analyse.info(system, integrator, per_atom=True)

# Run system with capped potentials, thermostat and increasing LJ epsilon
for k in range(args.warmup_loops):
    LJCapped.setPotential(0,0,espresso.interaction.LennardJonesCapped(epsilon_start + (epsilon-epsilon_start)*k*1.0/(args.warmup_loops-1), chain_setup.sigma, chain_setup.rc, caprad=chain_setup.caprad_LJ))
    integrator.run(args.warmup_steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)

# Remove LJ Capped potential
system.removeInteraction(0)

# Add non-capped LJ potential
LJ = espresso.interaction.VerletListLennardJones(verletList)
LJ.setPotential(0, 0, espresso.interaction.LennardJones(epsilon, chain_setup.sigma, chain_setup.rc))
system.addInteraction(LJ)

# Run system with non-capped potentials, thermostat and fixed LJ epsilon
for k in range(args.warmup_loops):
    integrator.run(args.warmup_steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)

thermostat.disconnect()
resetter = espresso.analysis.TotalVelocity(system)
resetter.reset()

fpl = espresso.FixedPairList(system.storage)
potMirrorLennardJones = espresso.interaction.MirrorLennardJones(epsilon=1.0, sigma=1.0)
interMirrorLennardJones = espresso.interaction.FixedPairListMirrorLennardJones(system, fpl, potMirrorLennardJones)
system.addInteraction(interMirrorLennardJones)

AR = espresso.integrator.AssociationReaction(system, verletList, fpl, system.storage)
AR.typeA = 0
AR.typeB = 0
AR.deltaA = -1
AR.deltaB = 2
AR.stateAMin = 2
AR.interval = args.interval
AR.rate = args.rate
AR.cutoff = chain_setup.rc
integrator.addExtension(AR)

integrator.step=0

import chain_h5md

if args.file is not None:
    traj_file = chain_h5md.DumpH5MD(args.file, system, integrator, 'Pierre de Buyl', edges=[edge_i for edge_i in system.bc.boxL], n_states=3)
    if args.stop_at<0:
        topo_file = chain_h5md.DumpTopo(traj_file, 'atoms', 'chains', system, integrator, fpl, time=True, chunks=(8,128,2))

last_dump_step=0
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon and crosslinking
for k in range(1,args.loops+1):
    integrator.run(args.steps)
    fpls = fpl.size()
    print fpls,
    espresso.tools.analyse.info(system, integrator, per_atom=True)
    if args.file is not None: traj_file.analyse()
    if args.file is not None and args.stop_at<0 and k%args.dump_interval==0:
        traj_file.dump()
        if args.stop_at<0: topo_file.dump()
    if args.stop_at>=0 and sum(fpls)>=args.stop_at:
        print "Chain size: ", sum(fpls)
        print "Last crosslinking step", integrator.step
        break

if args.stop_at<0:
    if args.file is not None: traj_file.close()
    import sys
    sys.exit()

# Disconnect AR
AR.disconnect()

# Store topology
if args.file is not None:
    chain_h5md.DumpTopo(traj_file, 'atoms', 'chains', system, integrator, fpl)

# Run system with non-capped potentials, no thermostat, fixed LJ epsilon and no
# crosslinking
for k in range(1, args.loops+1):
    integrator.run(args.steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)
    if args.file is not None: traj_file.analyse()
    if args.file is not None and k%args.dump_interval==0:
        traj_file.dump()

