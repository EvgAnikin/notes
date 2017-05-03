import sys
import math
import numpy as np
import random
import scipy.linalg as linalg
import matplotlib.pyplot as plt
from collections import namedtuple


def hamiltonian(N, t, epsilon):
    ham = np.zeros(N*N, dtype=complex).reshape((N,N))
    potential = 2*epsilon*(np.random.rand(N) - 0.5)
    for i in xrange(N-1):
        ham[i,i+1] = ham[i+1,i] = -t
    for i in xrange(N):
        ham[i,i] = potential[i]
    return ham


if __name__ == '__main__':
    energies, wfs = linalg.eigh(hamiltonian(N=400, t=1, epsilon=0.5))
    
    while True:
        print 'Enter the number of wavefunction'
        c = raw_input()
        if c == 'stop':
            break
        else:
            try:
                plt.plot(abs(wfs[:,int(c)]))
                plt.show()
            except IndexError:
                print 'Index out of range'
            except ValueError:
                print 'Can\'t read the number'
