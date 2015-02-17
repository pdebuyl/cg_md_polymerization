import espresso
import mpi4py.MPI as MPI

def get_velocity(system, n):
    """Obtain total velocity of a espresso system."""
    total_v = espresso.Real3D(0.)
    total_m = 0.
    for i in range(n):
        p = system.storage.getParticle(i)
        total_v += p.v*p.mass
        total_m += p.mass
    return total_v/total_m

def reset_velocity(system, n):
    """Reset the total velocity of a espresso system."""
    excess_v = get_velocity(system, n)
    for i in range(n):
        v = system.storage.getParticle(i).v
        system.storage.modifyParticle(i, 'v', v-excess_v)

# LJ settins
sigma = 1.0
epsilon=1.0
caprad_LJ=0.85
rc = pow(2., 1./6.)

# General settings
skin = 0.3

def monomer_system(num_particles, density=0.8, seed=None):

    num_particles = num_particles
    L = pow(num_particles/density, 1./3.)
    box = (L, L, L)

    # Initialize the espresso system
    system         = espresso.System()
    if seed is not None:
        system.rng     = espresso.esutil.RNG(seed)
    else:
        system.rng     = espresso.esutil.RNG()
    system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
    system.skin    = skin
    nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    def normal_v():
        return espresso.Real3D(system.rng.normal()*0.5, system.rng.normal()*0.5, system.rng.normal()*0.5)

    # Add the individual particles
    Xs = []
    pid = 0
    for i in range(num_particles):
        pos = system.bc.getRandomPos()
        v = espresso.Real3D(system.rng.normal(),system.rng.normal(),system.rng.normal())
        Xs.append([pid, pos, v])
        pid += 1
    system.storage.addParticles(Xs, 'id', 'pos', 'v')

    # Define capped LJ potential
    verletList = espresso.VerletList(system, cutoff=rc)
    LJCapped    = espresso.interaction.VerletListLennardJonesCapped(verletList)
    LJCapped.setPotential(type1=0, type2=0, potential=espresso.interaction.LennardJonesCapped(epsilon=epsilon, sigma=sigma, cutoff=rc, caprad=caprad_LJ))
    system.addInteraction(LJCapped)

    # Define integrator and StochasticVelocityRescaling thermostat
    integrator     = espresso.integrator.VelocityVerlet(system)
    thermostat = espresso.integrator.StochasticVelocityRescaling(system)
    thermostat.temperature = 1.0
    integrator.addExtension(thermostat)

    system.storage.decompose()

    return system, integrator, LJCapped, verletList, thermostat, num_particles
