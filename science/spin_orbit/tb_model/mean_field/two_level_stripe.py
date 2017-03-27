from collections import namedtuple
from itertools import izip, repeat
import math
from matplotlib import pyplot as plt
import numpy as np
from scipy.linalg import eigh, eigvalsh
#from tblib import *


def two_level_stripe_hamiltomian(py, args, N, spin, potential=None):
    xi, m, t = args
    H = np.zeros((2*N)**2, dtype=complex).reshape(N,2,N,2)

    if potential is None:
        potential = np.zeros(N)

    s = 1 if spin == 'up' else -1
    
    for i in xrange(N):
       H[i,:,i,:] = np.array([xi + 1./m*(2 - math.cos(py)) + potential[i], 
                              -2j*s*t*math.sin(py),
                              2j*s*t*math.sin(py), 
                              -xi - 1./m*(2 - math.cos(py)) 
                                   + potential[i]]).reshape(2,2)
       if i > 0:
            H[i-1,:,i,:] = np.array([0.5/m, -1j*t*s,
                                     -1j*t*s, -0.5/m]).reshape(2,2)
            H[i,:,i-1,:] = np.array([0.5/m, 1j*t*s,
                                     1j*t*s, -0.5/m]).reshape(2,2)
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


def get_wf_density(wf, band_index):
    return np.sum(np.absolute(wf)**2, axis=band_index)


def get_densities(ham_up, ham_down, nx, ny, n_particles):
    pyrange = np.linspace(0, 2*math.pi, ny, endpoint=False)

    State = namedtuple('State', ['energy', 'wf', 'py', 'spin']) 
    states = []
    for py in pyrange:
        energies_up, wfs_up = eigh(ham_up(py))
        energies_down, wfs_down = eigh(ham_down(py))

        states_up = [State(*s) for s in izip(energies_up, wfs_up.T, 
                                            repeat(py), repeat('up'))]
        states_down = [State(*s) for s in izip(energies_down, wfs_down.T, 
                                            repeat(py), repeat('down'))]

        states.extend(states_up)
        states.extend(states_down)

    states.sort(key=lambda state: state.energy)

    density_up = np.zeros(nx)
    density_down = np.zeros(nx)

    for s in states[:n_particles]:
        if s.spin == 'up': 
            density_up += get_wf_density(s.wf.reshape(nx, 2), band_index=1) 
        else:
            density_down += get_wf_density(s.wf.reshape(nx, 2), band_index=1)

    fermi_energy = states[n_particles-1].energy

    return density_up/ny, density_down/ny, fermi_energy


def update_densities(density_up, density_down, args, nx, ny, g, ex_pot, n_particles):
    ham_up = lambda py: two_level_stripe_hamiltomian(py, args, nx, 'up', 
                        g*density_down + ex_pot) 
    ham_down = lambda py: two_level_stripe_hamiltomian(py, args, nx, 'down', 
                        g*density_up + ex_pot) 
        
    d_up_new, d_down_new, fermi_energy = get_densities(ham_up, ham_down, nx, ny, n_particles)
    sm_f = 0.6
    return sm_f*density_up + (1 - sm_f)*d_up_new, \
           sm_f*density_down + (1 - sm_f)*d_down_new, fermi_energy


def get_external_potential(nx, length, strength):
    external_potential = np.zeros(nx)
    if length > 0:
        external_potential[:length] = strength*np.linspace(1,0,length)
        external_potential[-length:] = strength*np.linspace(0,1,length)
    return external_potential


def plot_energy_bands():
    energies = []
    pyrange = np.linspace(-math.pi, math.pi, 41)
    for py in pyrange:
        energies.append(
            eigvalsh(two_level_stripe_hamiltomian(py, (xi, m, t), nx, 'up', g*d_down)))
    energies = np.array(energies).T
    
    for band in energies:
        plt.plot(pyrange, band)
    plt.show()



if __name__ == '__main__':
    xi, m, t = -0.2, 1, 0.5 
    potential_length = 8
    potential_strength = 2
    nx, ny = 15 + 2*potential_length, 20
    n_iter = 10
    epsilon = 0.1
    g = 5
    n_part = int(round(2*nx*ny*0.97))

    d_up = np.ones(nx) + epsilon*np.random.rand(nx)
    d_down = np.ones(nx) + epsilon*np.random.rand(nx)
    ex_pot = get_external_potential(nx, potential_length, potential_strength)

    for i in xrange(n_iter):
        print '{}th iteration...'.format(i)
        d_up, d_down, fermi_energy = update_densities(d_up, d_down, (xi, m, t), 
                                                      nx, ny, g, ex_pot, n_part)

    density_figure = plt.figure(1)
    plt.plot(range(0,nx), d_up)
    plt.plot(range(0,nx), d_down)
    plt.xlim(0,nx-1)
#    plt.ylim(0.9,1.2)
    density_figure.show()

#=============================================================================================

    band_figure_up = plt.figure(2)


    energies = []
    Npy = 41
    pyrange = np.linspace(-math.pi, math.pi, Npy)
    for py in pyrange:
        energies.append(
            eigvalsh(two_level_stripe_hamiltomian(py, (xi, m, t), nx, 
                                                  'up', g*d_down + ex_pot)))
    energies = np.array(energies).T
    
    for band in energies:
        plt.plot(pyrange, band, color='black')

    plt.plot(pyrange, fermi_energy*np.ones(len(pyrange)), color='red')
    
    band_figure_up.show()

#=============================================================================================

    band_figure_down = plt.figure(3)
    energies = []
    pyrange = np.linspace(-math.pi, math.pi, 41)
    for py in pyrange:
        energies.append(
            eigvalsh(two_level_stripe_hamiltomian(py, (xi, m, t), nx, 
                                                  'down', g*d_up + ex_pot)))
    energies = np.array(energies).T
    
    for band in energies:
        plt.plot(pyrange, band, color='black')

    plt.plot(pyrange, fermi_energy*np.ones(len(pyrange)), color='red')

    band_figure_down.show()

    raw_input()
