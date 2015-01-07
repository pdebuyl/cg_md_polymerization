from __future__ import print_function, division

import numpy as np

N = 2001
eps = 0.1

rc = 2.**(1./6.)
r = np.linspace(0, 2*rc, N+2)[1:-1]

support = (r>=rc)

def LJ(r):
    rm6 = r**(-6)
    return 4*rm6*(rm6-1)+1

def F_LJ(r):
    return 24*(2*r**(-13)-r**(-7))

def FPRIME_LJ(r):
    return -24*(26*r**(-14)-7*r**(-8))

def mirror_LJ(r):
    rm6 = (2*rc - r)**(-6)
    return 4*rm6*(rm6-1)+1

def mirror_F_LJ(r):
    rr = 2*rc - r
    return -24*(2*rr**(-13)-rr**(-7))

def mirror_FPRIME_LJ(r):
    rr = 2*rc - r
    return 24*(26*rr**(-14)-7*rr**(-8))

V = np.zeros_like(r)
V[support] = mirror_LJ(r[support])

do_plot = False
if do_plot:
    import matplotlib.pyplot as plt

    plt.plot(r, LJ(r))
    plt.plot(r, F_LJ(r), ls='--')
    plt.plot(r, V)
    plt.plot(r[support], mirror_F_LJ(r[support]), ls='--')
    plt.ylim(-2,5)
    plt.show()

r = np.linspace(eps, 2*rc-eps, N)


print("# rlo={rlo}, rhi={rhi}".format(rlo=rc, rhi=2*rc-eps))
print("""# Mirror Lennard-Jones potential

MIRROR_LJ
N {N} FP {fplo} {fphi}
""".format(N=N, rlo=eps, rhi=2*rc-eps, fplo=FPRIME_LJ(eps), fphi=mirror_FPRIME_LJ(2*rc-eps)))

for idx, x in enumerate(r):
    if x<rc:
        print("%i %f %f %f" % (idx+1, x, LJ(x), F_LJ(x)))
    else:
        print("%i %f %f %f" % (idx+1, x, mirror_LJ(x), mirror_F_LJ(x)))
