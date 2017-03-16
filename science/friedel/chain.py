import sys
import math
import numpy as np
import random
import scipy.linalg as linalg
import matplotlib.pyplot as plt
from collections import namedtuple


def chain_hamiltonian(N, t, g, density=None):
    ham = np.zeros(N*N).reshape((N,N))
    for i in xrange(N-1):
        ham[i,i+1] = ham[i+1,i] = -t
    if density is not None:
        for i in xrange(N):
            ham[i,i] = g*density[i]
    return ham


def get_densities(ham_up, ham_down, n_particles):
    energies_up, wf_up = linalg.eigh(ham_up)
    energies_down, wf_down = linalg.eigh(ham_down)


    N = ham_up.shape[0]
    State = namedtuple('State', ['energy', 'wf', 'spin']) 
    states_up = [State(*s) for s in zip(energies_up, wf_up.T, ['up']*N)]
    states_down = [State(*s) for s in zip(energies_down, wf_down.T, ['down']*N)]

    states = states_up + states_down
    states.sort(key = lambda state: state.energy)

    density_up = np.zeros(N)
    density_down = np.zeros(N)
    for i in xrange(n_particles):
        state = states[i]
        if state.spin == 'up':  
            density_up += state.wf*state.wf
        else:
            density_down += state.wf*state.wf
    
    return density_up, density_down


def update_densities(density_up, density_down, N, t, g, n_particles):
    tau = 0.5
    density_up_new, density_down_new = get_densities(
                    chain_hamiltonian(N, t, g, density_down), 
                    chain_hamiltonian(N, t, g, density_up),
                    n_particles)
    return (1 - tau)*density_up + tau*density_up_new,\
           (1 - tau)*density_down + tau*density_down_new


if __name__ == '__main__':
    N = 100
    t = 1
    g = 1
    n_part = 20
    n_iter = 100
    epsilon = 0.01

#    stoner_parameter = g/(2*math.pi*t*math.sin(math.pi*n_part/N))
#    print stoner_parameter

    random.seed()
#    d_up = epsilon*np.random.rand(N)
#    d_down = epsilon*np.random.rand(N)
    d_up = np.concatenate((epsilon*np.ones(N/2), -epsilon*np.ones(N/2)))
    d_down = np.concatenate((-epsilon*np.ones(N/2), epsilon*np.ones(N/2)))

    for i in range(n_iter):
        d_up, d_down = update_densities(d_up, d_down, N, t, g, n_part)

    plt.plot(range(N), d_up)
    plt.plot(range(N), d_down)
    plt.xlim(-1,N)
    plt.show()
