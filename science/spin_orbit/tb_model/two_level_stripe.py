import math
from matplotlib import pyplot as plt
import numpy as np
from scipy.linalg import eigh
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


def plot_edge_state(py, number, N, xi, m, t):
    state_figure = plt.figure(1)
    ham = two_level_stripe_hamiltomian(py=py, N=N, xi=xi, m=m, t=t)
    energies, states_raw = eigh(ham)
    states = states_raw.reshape(N,2,2*N).transpose(2,0,1)
    densities = np.sum(abs(states)**2, axis=2)
    plt.plot(densities[number])
    state_figure.show()


def plot_bands(N, xi, m, t):
    plim = (-math.pi,math.pi)
    pyrange, energies = stripe_energies(two_level_stripe_hamiltomian, 
                                       (N, xi, m, t),
                                       plim = plim,
                                       NX = 51)
    band_figure = plt.figure(0)
    for level in energies:
        plt.plot(pyrange,level)
    plt.xlim(*plim)
    band_figure.show()


if __name__ == '__main__':
    N = 100
    xi, m, t = -0.2, 1, 0.4
    
    plot_bands(N, xi, m, t)
#    plot_edge_state(0.0, N, N, xi, m, t)
