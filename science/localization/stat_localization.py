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


def mean_x(wf):
    sum_x = 0
    for i in xrange(len(wf)):
        sum_x += i*abs(wf[i])**2
    return sum_x


def disp_x(wf):
    x = mean_x(wf)
    sum_sq_dev = 0
    for i in xrange(len(wf)):
       sum_sq_dev += (i - x)**2*abs(wf[i])**2
    return math.sqrt(sum_sq_dev)


if __name__ == '__main__':
    energies, wfs = linalg.eigh(hamiltonian(N=800, t=1, epsilon=0.5))
    wfs = wfs.T
    dispersions = [disp_x(wf) for wf in wfs]

    fig1 = plt.figure(1)
    plt.plot(energies, dispersions)
    plt.show()
    
<<<<<<< HEAD
#    while True:
#        print 'Enter the number of wavefunction'
#        c = raw_input()
#        if c == 'stop':
#            break
#        else:
#            try:
#                plt.plot(wfs[int(c)])
#                plt.show()
#            except IndexError:
#                print 'Index out of range'
#            except ValueError:
#                print 'Can\'t read the number'
=======
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
>>>>>>> 3d4448b3e276ff1b891a38c21badedc55c202bfd
