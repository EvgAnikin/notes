import math
from matplotlib import pyplot as plt
import numpy as np
from tblib import *


def two_level_stripe_hamiltomian(py, N, xi, m, t):
    H = np.zeros((2*N)**2, dtype=complex).reshape(N,2,N,2)
    
    voltage = 0.2
    potential = np.zeros(N, dtype=complex)

#    for i in xrange(N):
#       potential[i] = -voltage*(i - N/2.)/N
    xi_local = xi
    for i in xrange(N):
       H[i,:,i,:] = np.array([xi_local + 1./m*(2 - math.cos(py)) + potential[i], 
                              -2j*t*math.sin(py),
                              2j*t*math.sin(py), 
                              -xi_local - 1./m*(2 - math.cos(py)) 
                                   + potential[i]]).reshape(2,2)
       if i > 0:
            H[i-1,:,i,:] = np.array([0.5/m, -1j*t,
                                     -1j*t, -0.5/m]).reshape(2,2)
            H[i,:,i-1,:] = np.array([0.5/m, 1j*t,
                                     1j*t, -0.5/m]).reshape(2,2)
    return H.reshape(2*N,2*N)


if __name__ == '__main__':
    plim = (-math.pi,math.pi)
    pyrange, energies = stripe_energies(two_level_stripe_hamiltomian, 
                                       (50, -0.1, 1, 0.5),
                                       plim = plim,
                                       NX = 101)
    for level in energies:
        plt.plot(pyrange,level)
    plt.xlim(*plim)
#    plt.ylim(-1,1)
    plt.show()
    
