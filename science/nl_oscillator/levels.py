#encoding:utf-8
from __future__ import print_function
import math
import numpy as np
from numpy.linalg import eig, eigh, eigvalsh, lstsq
import matplotlib.pyplot as plt


def nonlinear_oscillator(omega, beta, E, N):
    ham = np.zeros(N**2).reshape(N,N)

    for i in xrange(N):
        ham[i,i] = omega*i + beta/2*(i**2 - 1)
        if i < N - 1:
            ham[i,i+1] = math.sqrt(i+1)*E
            ham[i+1,i] = math.sqrt(i+1)*E

    return ham


N_E = 400
omega_0=1
Omega=1.025
beta=0.003
N = 60
ext_fields = np.linspace(0, 0.2, N_E)

for E in ext_fields:
    energies = eigvalsh(nonlinear_oscillator(omega_0 - Omega, beta, E, N))
    plt.scatter([E]*len(energies), energies, marker='.', s=1)

plt.xlim(0, 0.05)
plt.ylim(-1,1)
plt.show()
