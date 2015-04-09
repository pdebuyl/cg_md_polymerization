import argparse

parser = argparse.ArgumentParser()
parser.add_argument('nchains', type=int, help='Number of polymer chains')
parser.add_argument('nx', type=int, help='Number of crosslinker units')
parser.add_argument('--functionality', type=int, help='Functionality of crosslinker.', default=5)
parser.add_argument('--chain-length', type=int, help='Chain length.', default=5)
parser.add_argument('--dt', type=float, help='Timestep.', default=0.01)
parser.add_argument('--warmup-steps', type=int, help='Number of steps in warmup loop.', default=100)
parser.add_argument('--warmup-loops', type=int, help='Number of iterations in warmup loop.', default=100)
parser.add_argument('--steps', type=int, help='Number of steps in main loop.', default=100)
parser.add_argument('--loops', type=int, help='Number of iterations in main loop.', default=100)
parser.add_argument('--post-loops', type=int, help='Number of iterations in post loop.', default=-1)
parser.add_argument('--density', type=float, help='Number density.', default=0.8)
parser.add_argument('--rate', type=float, help='Reaction rate.', default=0.1)
parser.add_argument('--interval', type=int, help='Steps between reactions.', default=100)
parser.add_argument('--file', type=str, help='Filename for H5MD output.', default=None)
parser.add_argument('--dump-interval', type=int, help='Interval at which to dump the configuration.', default=20)
parser.add_argument('--seed', type=int, help='Seed for the rng.')
args = parser.parse_args()

import espressopp
import epoxy_setup

system, integrator, LJCapped, verletList, FENECapped, chainFPL, thermostat, num_particles = epoxy_setup.chains_x_system(args.nchains, args.chain_length, args.nx, density=args.density, seed=args.seed)

integrator.dt = args.dt

epsilon_start = 0.1
epsilon = 1.0

espressopp.tools.analyse.info(system, integrator, per_atom=True)

# Run system with capped potentials, thermostat and increasing LJ epsilon
for k in range(args.warmup_loops):
    LJCapped.setPotential(0,0,espressopp.interaction.LennardJonesCapped(epsilon_start + (epsilon-epsilon_start)*k*1.0/(args.warmup_loops-1), epoxy_setup.sigma, epoxy_setup.rc, caprad=epoxy_setup.caprad_LJ))
    integrator.run(args.warmup_steps)
    espressopp.tools.analyse.info(system, integrator, per_atom=True)

# Remove FENE Capped potential
system.removeInteraction(1)
# Remove LJ Capped potential
system.removeInteraction(0)

# Add non-capped LJ potential
LJ = espressopp.interaction.VerletListLennardJones(verletList)
LJ.setPotential(0, 0, espressopp.interaction.LennardJones(epsilon, epoxy_setup.sigma, epoxy_setup.rc))
system.addInteraction(LJ)

# Add non-capped FENE potential
FENE = espressopp.interaction.FixedPairListFENE(system, chainFPL, espressopp.interaction.FENE(epoxy_setup.K, 0.0, epoxy_setup.rMax))
system.addInteraction(FENE)

# Run system with non-capped potentials, thermostat and fixed LJ epsilon
for k in range(args.warmup_loops):
    integrator.run(args.warmup_steps)
    espressopp.tools.analyse.info(system, integrator, per_atom=True)

thermostat.disconnect()
resetter = espressopp.analysis.TotalVelocity(system)
resetter.reset()

for i in range(args.nchains):
  for j in range(1,args.chain_length-1):
    system.storage.modifyParticle(i*args.chain_length+j, 'state', 1)
pid=args.nchains*args.chain_length
for i in range(args.nx):
  system.storage.modifyParticle(pid, 'type', 1)
  system.storage.modifyParticle(pid, 'state', args.functionality)
  pid += 1

LJ.setPotential(0, 1, espressopp.interaction.LennardJones(epsilon, epoxy_setup.sigma, epoxy_setup.rc))
LJ.setPotential(1, 1, espressopp.interaction.LennardJones(epsilon, epoxy_setup.sigma, epoxy_setup.rc))

fpl = espressopp.FixedPairList(system.storage)
potMirrorLennardJones = espressopp.interaction.MirrorLennardJones(epsilon=1.0, sigma=1.0)
interMirrorLennardJones = espressopp.interaction.FixedPairListMirrorLennardJones(system, fpl, potMirrorLennardJones)
system.addInteraction(interMirrorLennardJones)

AR = espressopp.integrator.AssociationReaction(system, verletList, fpl, system.storage)
AR.typeA = 1
AR.typeB = 0
AR.deltaA = -1
AR.deltaB = 1
AR.stateAMin = 1
AR.interval = args.interval
AR.rate = args.rate
AR.cutoff = epoxy_setup.rc
integrator.addExtension(AR)

integrator.step=0

import epoxy_h5md

if args.file is not None:
    traj_file = epoxy_h5md.DumpH5MD(args.file, system, integrator, 'Pierre de Buyl', edges=[edge_i for edge_i in system.bc.boxL], n_states=args.functionality+1, species_to_count=1)
    epoxy_h5md.DumpTopo(traj_file, 'atoms', 'chains', system, integrator, chainFPL)
    if args.post_loops<0:
        topo_file = epoxy_h5md.DumpTopo(traj_file, 'atoms', 'crosslinks', system, integrator, fpl, time=True, chunks=(8,128,2))

# Run system with non-capped potentials, no thermostat, fixed LJ epsilon and crosslinking
if args.file is not None: traj_file.analyse()
for k in range(1, args.loops+1):
    integrator.run(args.steps)
    print fpl.size(),
    espressopp.tools.analyse.info(system, integrator, per_atom=True)
    if args.file is not None and k%2==0: traj_file.analyse()
    if args.file is not None and args.post_loops<0 and k%args.dump_interval==0:
        traj_file.dump()
        topo_file.dump()

if args.post_loops<0:
    if args.file is not None: traj_file.close()
    import sys
    sys.exit()

# Disconnect AR
AR.disconnect()
if args.file is not None: epoxy_h5md.DumpTopo(traj_file, 'atoms', 'crosslinks', system, integrator, fpl)

# Run system with non-capped potentials, no thermostat, fixed LJ epsilon and crosslinking
for k in range(1, args.post_loops+1):
    integrator.run(args.steps)
    print fpl.size(),
    espressopp.tools.analyse.info(system, integrator, per_atom=True)
    if args.file is not None and k%args.dump_interval==0:
        traj_file.dump()

if args.file is not None: traj_file.close()
