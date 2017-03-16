from collections import namedtuple
import math
from matplotlib import pyplot as plt
import numpy as np
from scipy.linalg import eigh
#from tblib import *


def two_level_stripe_hamiltomian(py, args, spin, density=None):
    N, xi, m, t, g = args
    H = np.zeros((2*N)**2, dtype=complex).reshape(N,2,N,2)

    if density is None:
        density = np.zeros(N)

    s = 1 if spin == 'up' else -1
    
    for i in xrange(N):
       H[i,:,i,:] = np.array([xi + 1./m*(2 - math.cos(py)) + density[i], 
                              -2j*s*t*math.sin(py),
                              2j*s*t*math.sin(py), 
                              -xi - 1./m*(2 - math.cos(py)) 
                                   + density[i]]).reshape(2,2)
       if i > 0:
            H[i-1,:,i,:] = np.array([0.5/m, -1j*t,
                                     -1j*t, -0.5/m]).reshape(2,2)
            H[i,:,i-1,:] = np.array([0.5/m, 1j*t,
                                     1j*t, -0.5/m]).reshape(2,2)
    return H.reshape(2*N,2*N)


def primitive_stripe(py, t, N, potential=None):
    H = np.zeros(N**2, dtype=complex).reshape(N,N)

    if potential is None:
        potential = np.zeros(N) 

    for i in xrange(N):
        H[i,i] = -2*t*math.cos(py) + potential[i]
        if i < N - 1:
            H[i,i+1] = H[i+1,i] = -t

    return H


def get_densities(ham_up, ham_down, nx, ny):
    pyrange = np.linspace(0, 2*math.pi, ny, endpoint=False)

    State = namedtuple('State', ['energy', 'wf', 'py', 'spin']) 
    states = []
    for py in pyrange:
        energies_up, wfs_up = eigh(ham_up(py))
        energies_down, wfs_down = eigh(ham_down(py))

        states_up = [State(*s) for s in zip(energies_up, wfs_up.T, [py]*nx, ['up']*nx)]
        states_down = [State(*s) for s in zip(energies_down, wfs_down.T, [py]*nx, ['down']*nx)]

        states.extend(states_up)
        states.extend(states_down)

    states.sort(key=lambda state: state.energy)

    density_up = np.zeros(nx)
    density_down = np.zeros(nx)

    average_density = 0.018
    n_particles = 2*int(round(nx*ny*average_density)) # !
    
    for s in states[:n_particles]:
        if s.spin == 'up':
            density_up += np.absolute(s.wf)**2
        else:
            density_down += np.absolute(s.wf)**2

    return density_up/ny, density_down/ny


#def update_densities(density_up, density_down, t, g, N):
#    h_up = lambda py: primitive_stripe(py, t, N, potential=g*density_down)
        


if __name__ == '__main__':
    t = 1
    nx = 50
    ny = 19
    ham_up = ham_down = lambda py: primitive_stripe(py, t, nx)

    d_up, d_down = get_densities(ham_up, ham_down, nx, ny)

    plt.plot(range(1,nx+1), d_up)
    plt.show()
