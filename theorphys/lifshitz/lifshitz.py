from __future__ import print_function
import math
import numpy as np
from numpy.linalg import eigvalsh
import matplotlib.pyplot as plt


def hamiltonian(n, L, a):
    N = int(n*L**3)
    nodes = []
    for i in xrange(N):
        x = np.random.rand()*L
        y = np.random.rand()*L
        z = np.random.rand()*L
        nodes.append((x,y,z))

    H = np.zeros(N*N).reshape(N, N)

    for i in xrange(N):
        for j in xrange(N):
            x1, y1, z1 = nodes[i]
            x2, y2, z2 = nodes[j]
            r_ij = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            if i != j:
                H[i,j] = math.exp(-r_ij/a)

    return H


energies = []
for i in xrange(10):
    print(i)
    H = hamiltonian(0.003, 50, 1)
    energies.append(eigvalsh(H))

hist, bins = np.histogram(energies, bins=100)

plt.plot(0.5*(bins[:-1] + bins[1:]), hist)
plt.show()
