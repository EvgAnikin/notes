import fractions
from itertools import product
import matplotlib.pyplot as plt
import math
import numpy as np
from numpy.linalg import eigvalsh


def hofstadter_ham(px, py, p, q):
    ham = np.zeros(q**2, dtype=complex).reshape((q,q))
    for i in xrange(q):
        ham[i,i] = 2*math.cos(py - i*2*math.pi*p/q)
        if i > 0:
            ham[i-1,i] = 1
            ham[i,i-1] = 1
    ham[q-1,0] = np.exp(1j*px)
    ham[0,q-1] = np.exp(-1j*px)
    return ham


if __name__ == '__main__':
    NX = NY = 4
    pxgrid = np.linspace(0,2*math.pi, NX)
    pygrid = np.linspace(0,2*math.pi, NY)

    q = 1000
    for p in xrange(q):
        f = fractions.Fraction(p,q)
        print 'diagonalising at value {}'.format(f)
        for px in pxgrid:
            for py in pygrid:
                H = hofstadter_ham(px,py,f.numerator, f.denominator)
                energies = eigvalsh(H)
                plt.scatter(len(energies)*[float(p)/q], energies, s = 0.1)
    plt.xlim(0,1)
    plt.show()
