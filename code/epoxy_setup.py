import espressopp
import mpi4py.MPI as MPI

def get_velocity(system, n):
    """Obtain total velocity of a espressopp system."""
    total_v = espressopp.Real3D(0.)
    total_m = 0.
    for i in range(n):
        p = system.storage.getParticle(i)
        total_v += p.v*p.mass
        total_m += p.mass
    return total_v/total_m

def reset_velocity(system, n):
    """Reset the total velocity of a espressopp system."""
    excess_v = get_velocity(system, n)
    for i in range(n):
        v = system.storage.getParticle(i).v
        system.storage.modifyParticle(i, 'v', v-excess_v)

# LJ settins
sigma = 1.0
epsilon=1.0
caprad_LJ=0.85
rc = pow(2., 1./6.)

# FENE settings
K=30.
rMax=1.5
caprad_FENE=1.4

# Polymer chain settings
bondlen=0.97

# General settings
skin = 0.3

def chains_x_system(num_chains, monomers_per_chain, num_X, density=0.8, seed=None):

    num_particles = num_chains*monomers_per_chain + num_X
    L = pow(num_particles/density, 1./3.)
    box = (L, L, L)

    # Initialize the espressopp system
    system         = espressopp.System()
    if seed is not None:
        system.rng     = espressopp.esutil.RNG(seed)
    else:
        system.rng     = espressopp.esutil.RNG()
    system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin    = skin
    nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    def normal_v():
        return espressopp.Real3D(system.rng.normal()*0.5, system.rng.normal()*0.5, system.rng.normal()*0.5)

    # Add the chains
    chainFPL = espressopp.FixedPairList(system.storage)
    pid      = 0
    for i in range(num_chains):
        chain=[]
        startpos = system.bc.getRandomPos()
        positions, bonds = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
        for k in range(monomers_per_chain):  
            part = [pid + k, positions[k], normal_v()]
            chain.append(part)
        pid += monomers_per_chain
        system.storage.addParticles(chain, 'id', 'pos', 'v')
        chainFPL.addBonds(bonds)

    # Add the individual particles
    Xs = []
    for i in range(num_X):
        pos = system.bc.getRandomPos()
        v = espressopp.Real3D(system.rng.normal(),system.rng.normal(),system.rng.normal())
        Xs.append([pid, pos, v])
        pid += 1
    system.storage.addParticles(Xs, 'id', 'pos', 'v')

    # Define capped LJ potential
    verletList = espressopp.VerletList(system, cutoff=rc)
    LJCapped    = espressopp.interaction.VerletListLennardJonesCapped(verletList)
    LJCapped.setPotential(type1=0, type2=0, potential=espressopp.interaction.LennardJonesCapped(epsilon=epsilon, sigma=sigma, cutoff=rc, caprad=caprad_LJ))
    system.addInteraction(LJCapped)

    # Define capped FENE potential
    potFENE   = espressopp.interaction.FENECapped(K=K, r0=0.0, rMax=rMax, caprad=caprad_FENE)
    FENECapped = espressopp.interaction.FixedPairListFENECapped(system, chainFPL, potFENE)
    system.addInteraction(FENECapped)

    # Define integrator and StochasticVelocityRescaling thermostat
    integrator     = espressopp.integrator.VelocityVerlet(system)
    thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
    thermostat.temperature = 1.0
    integrator.addExtension(thermostat)

    system.storage.decompose()

    return system, integrator, LJCapped, verletList, FENECapped, chainFPL, thermostat, num_particles
