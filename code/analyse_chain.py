#!/usr/bin/env python

import sys
import os
import os.path

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--lammps', type=str, nargs='+',
                    help='directories containing LAMMPS simulation files', default=[])
parser.add_argument('--espp', type=str, nargs='+',
                    help='directories containing ESPResSo++ simulation files', default=[])
parser.add_argument('--rate', type=float, default=0.1)
parser.add_argument('--sites', type=int, default=1)
parser.add_argument('-N', type=int, default=10000)
parser.add_argument('--fit', action='store_true')
parser.add_argument('--plot-all', action='store_true')
args = parser.parse_args()

import numpy as np
from scipy.optimize import leastsq
from io import StringIO
import matplotlib.pyplot as plt
import h5py

plt.rcParams['figure.figsize']= (10,6)
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 2

NNEIGH=3.25
fraction = float(args.sites)/float(args.N)

# Open lammps log file to extract thermodynamic observables
def from_log(logfile,i0,i1):
    return np.loadtxt(StringIO(u''.join(logfile[i0+1:i1])), unpack=True)

fitfunc = lambda p, t: 1*(1.-np.exp(-t*p[0]-p[1]))
errfunc = lambda p, t, y: fitfunc(p, t) - y

time_max = 0

p_data = []
phi_data = []
for d in args.lammps:
    logfile = open(os.path.join(os.getcwd(), d, 'log.lammps')).readlines()
    start_indices = [(i,l) for (i,l) in enumerate(logfile) if l.startswith('Time ')]
    stop_indices = [(i,l) for (i,l) in enumerate(logfile) if l.startswith('Loop time')]
    time, e_tot, temp, e_kin, e_vdw, e_bond, e_pot, press, rho, n_bonds, n_bonds_max, bonds = from_log(logfile, start_indices[-1][0], stop_indices[-1][0])
    time -= time[0]
    time_max = max(time_max, time.max())
    n_bonds += fraction
    phi_data.append(n_bonds)
    if args.plot_all: plt.plot(time, n_bonds)
    p, success = leastsq(errfunc, [args.rate*NNEIGH*fraction, 0./args.rate], args=(time, n_bonds))
    p_data.append(p)
    print p

if len(args.lammps)>0:
    plt.plot(time, np.mean(phi_data, axis=0))
    p_data = np.array(p_data)
    print "lmp fit rate", p_data.mean(axis=0)[0]
    if args.fit:
        plt.plot(time, fitfunc(p_data.mean(axis=0), time), 'k--')

p_data = []
phi_data = []
for d in args.espp:
    a = h5py.File(os.path.join(os.getcwd(), d, 'dump.h5'), 'r')
    sc_time = a['/observables/statecount/time'][:]
    time_max = max(time_max, sc_time.max())
    sc = a['/observables/statecount/value']
    phi = 1. - sc[:,0]/float(args.N)
    phi_data.append(phi)
    p, success = leastsq(errfunc, [args.rate*NNEIGH*fraction, 0./args.rate], args=(sc_time[:], phi))
    if args.plot_all: plt.plot(sc_time, phi)
    a.close()
    p_data.append(p)
    print p

if len(args.espp)>0:
    plt.plot(sc_time, np.mean(phi_data, axis=0))
    p_data = np.array(p_data)
    print "esp fit rate", p_data.mean(axis=0)[0]
    if args.fit:
        plt.plot(time, fitfunc(p_data.mean(axis=0), time), 'k--')

print "th. rate", args.rate*NNEIGH*fraction
if time_max>0:
    time = np.linspace(0,1,512)*time_max
    plt.plot(time, 1*(1.-np.exp(-time*args.rate*NNEIGH*fraction)))

plt.show()
