#!/usr/bin/env python

import os
import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dirs', type=str, nargs='+',
                    help='directories containing simulation files')
parser.add_argument('--rate', type=float, default=0.1)
parser.add_argument('--sites', type=int, default=1)
parser.add_argument('-N', type=int, default=10000)
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.optimize import leastsq

NNEIGH=3.5
fraction = float(args.sites)/float(args.N)

fitfunc = lambda p, t: 1*(1.-np.exp(-t*p[0]-p[1]))
errfunc = lambda p, t, y: fitfunc(p, t) - y

p_data = []
for d in args.dirs:
    a = h5py.File(os.path.join(os.getcwd(), d, 'dump.h5'), 'r')
    sc_time = a['/observables/statecount/time'][:]
    sc = a['/observables/statecount/value']
    p, success = leastsq(errfunc, [args.rate*NNEIGH*fraction, 0./args.rate], args=(sc_time[:], sc[:,1]/float(args.N)))
    plt.plot(sc_time, sc[:,1]/float(args.N))
    a.close()
    p_data.append(p)
    print p

plt.plot(sc_time, 1*(1.-np.exp(-sc_time*args.rate*NNEIGH*fraction)))
p_data = np.array(p_data)
print "fit rate", p_data.mean(axis=0)[0]
print "th. rate", args.rate*NNEIGH*fraction
plt.plot(sc_time, fitfunc(p_data.mean(axis=0), sc_time), 'k--')

plt.show()
