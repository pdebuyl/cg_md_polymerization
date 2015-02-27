# Pierre de Buyl 2014
# This file is licensed under the modified BSD license

import espresso
from espresso import Real3D
import pyh5md
import numpy as np

def DumpH5MD(filename, system, integrator, author, author_email=None, edges=None, edges_time=False, n_states=None, species_to_count=1):
    espresso.Version().info()
    f = pyh5md.H5MD_File(filename, 'w', creator='espressopp', creator_version=espresso.Version().info(), author=author, author_email=author_email)
    atoms = f.particles_group('atoms')
    maxParticleID = int(espresso.analysis.MaxPID(system).compute())
    pos = atoms.trajectory('position', (maxParticleID+1,3), np.float64)
    species = atoms.trajectory('species', (maxParticleID+1,), np.int32)
    state = atoms.trajectory('state', (maxParticleID+1,), np.int32)
    if edges_time:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], edges=edges, time=True)
    else:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], edges=edges)
    f.NPart = espresso.analysis.NPart(system).compute()
    f.observable('particle_number', data=int(f.NPart), time=False)
    f.f['observables'].attrs['dimension']=3
    obs_dict = {}
    f.n_states = n_states
    f.species_to_count = species_to_count
    state_tuple = (('statecount', (n_states,), np.int32),) if n_states is not None else ()
    for o in (
            ('temperature', (), np.float64), ('kinetic_energy', (), np.float64),
            ('pressure', (), np.float64), ('pressure_tensor', (6,), np.float64),
            ('potential_energy', (), np.float64),('internal_energy', (), np.float64),
            ('lennard_jones', (), np.float64), ('fene', (), np.float64),
            ('mirror_lennard_jones', (), np.float64),
    ) + state_tuple:
        obs_dict[o[0]] = f.observable(*o)
    def dump():
        step = integrator.step
        time = integrator.step*integrator.dt
        maxParticleID = int(espresso.analysis.MaxPID(system).compute())
        if maxParticleID>pos.value.shape[1]:
            raise ValueError('System too large for dataset')
        r = np.array(
            [[x for x in system.storage.getParticle(pid).pos] for pid in range(maxParticleID+1)]
            )
        pos.append(r, step, time)
        del r
        species.append(
            np.array([system.storage.getParticle(pid).type for pid in range(maxParticleID+1)]),
            step, time)
        state.append(
            np.array([system.storage.getParticle(pid).state for pid in range(maxParticleID+1)]),
            step, time)
    f.dump = dump
    def analyse():
        step = integrator.step
        time = integrator.step*integrator.dt
        T = espresso.analysis.Temperature(system).compute()
        obs_dict['temperature'].append(T, step, time)
        P      = espresso.analysis.Pressure(system).compute()
        obs_dict['pressure'].append(P, step, time)
        Pij    = espresso.analysis.PressureTensor(system).compute()
        obs_dict['pressure_tensor'].append(Pij, step, time)
        Ek     = (3.0/2.0) * T
        obs_dict['kinetic_energy'].append(Ek, step, time)
        LJ = system.getInteraction(0).computeEnergy()/f.NPart
        obs_dict['lennard_jones'].append(LJ, step, time)
        FENE = system.getInteraction(1).computeEnergy()/f.NPart
        obs_dict['fene'].append(FENE, step, time)
        MirrorLJ = system.getInteraction(2).computeEnergy()/f.NPart
        obs_dict['mirror_lennard_jones'].append(MirrorLJ, step, time)
        potential = LJ+FENE+MirrorLJ
        obs_dict['potential_energy'].append(potential, step, time)
        obs_dict['internal_energy'].append(Ek+potential, step, time)
        if f.n_states is not None:
            local_state = np.array([system.storage.getParticle(pid).state for pid in range(maxParticleID+1) if system.storage.getParticle(pid).type==f.species_to_count])
            obs_dict['statecount'].append(np.bincount(local_state, minlength=f.n_states), step, time)
    f.analyse = analyse
    return f

class DumpTopo(object):
    def __init__(self, f, group, name, system, integrator, fpl, time=False, chunks=None):
        self.file = f
        self.system = system
        self.integrator = integrator
        self.fpl = fpl
        f.f.require_group('topology')
        f.f.require_group('topology/'+group)
        if time:
            self.element = pyh5md.base.TimeData(f.f['topology/atoms'], name, shape=(0,2), dtype=np.int32, chunks=chunks, compression='gzip', fillvalue=-1)
        else:
            data = np.array([b for local_bonds in fpl.getBonds() for b in local_bonds])
            self.element = pyh5md.base.FixedData(f.f['topology/atoms'], name, data=data)

    def dump(self):
        if not isinstance(self.element, pyh5md.base.TimeData):
            raise UserWarning("Trying to append data to a non suitable object.")
        bl = np.array([b for local_bonds in self.fpl.getBonds() for b in local_bonds])
        self.element.append(bl, self.integrator.step, self.integrator.step*self.integrator.dt)

