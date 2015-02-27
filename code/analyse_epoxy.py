#!/usr/bin/env python

import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--lammps', type=str, nargs='+',
                    help='directories containing LAMMPS simulation files', default=[])
parser.add_argument('--espp', type=str, nargs='+',
                    help='directories containing ESPResSo++ simulation files', default=[])
parser.add_argument('--rate', type=float, default=0.1)
parser.add_argument('--NC', type=int, default=2500)
parser.add_argument('--NX', type=int, default=1000)
parser.add_argument('--dt', type=float, default=0.0025)
parser.add_argument('--dump-interval', type=int, default=200)
args = parser.parse_args()

from io import StringIO
import gzip
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rcParams['figure.figsize']= (10,6)
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 2

def get_file(f):
    """Iterator over the timeframe in a lammps dump"""
    sio = StringIO()
    first = True
    finish = False
    while True:
        l = f.readline()
        if len(l)==0:
            finish=True
        if l.strip()=='ITEM: TIMESTEP' or finish:
            if first:
                first = False
            else:
                sio.seek(0)
                yield sio
                if finish: return
                sio = StringIO()
        sio.write(unicode(l))

for d in args.lammps:
    zf = gzip.open(os.path.join(d,'nb.txt.gz'), 'r')
    nb = []
    for onef in get_file(zf):
        nb.append(np.bincount(np.loadtxt(onef, skiprows=9, unpack=True, dtype=int), minlength=6))
    nb = np.array(nb)/float(args.NX)
    k_time = np.arange(nb.shape[0])*(args.dt*args.dump_interval*args.rate)
    plt.plot(k_time,nb)
    zf.close()

for d in args.espp:
    a = h5py.File(os.path.join(d, 'dump.h5'), 'r')
    sc_time = a['/observables/statecount/time'][:]
    sc = a['/observables/statecount/value'][:]/float(args.NX)
    plt.plot(sc_time*args.rate, sc)
    a.close()


plt.xlabel(r'$k t$')
plt.show()
