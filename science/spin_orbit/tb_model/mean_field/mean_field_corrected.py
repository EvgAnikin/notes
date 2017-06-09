from collections import namedtuple
from itertools import izip, repeat
import math
from matplotlib import pyplot as plt
import numpy as np
from scipy.linalg import eigh, eigvalsh
#from tblib import *


def get_external_potential(nx, length, strength):
    external_potential = np.zeros(nx)
    if length > 0:
        external_potential[:length] = strength*np.linspace(1,0,length)
        external_potential[-length:] = strength*np.linspace(0,1,length)
    return external_potential


def get_external_xi(nx, abs_xi, length):
    middle = divmod(nx,2)[0]
    ext_xi = np.zeros(nx)
    ext_xi[:middle - length] = np.ones(middle - length)*abs_xi
    for i in xrange(0, 2*length):
        ext_xi[middle - length + i] = abs_xi - 0.5*i/length
    ext_xi[middle + length:] = -np.ones(nx - middle - length)*abs_xi
    return ext_xi


class TwoLevelStripeHam:
    def get_hartree_fock_potential(self, spin): 
        opposite_spin = 'up' if spin == 'down' else 'down'
        return self.ext_potential + self.g_pot*self.densities[opposite_spin]


#    def get_hartree_fock_xi(self, spin):
#        opposite_spin = 'up' if spin == 'down' else 'down'
#        return self.ext_xi + self.g_xi*(self.density[opposite_spin] - np.ones(self.nx))


    def __init__(self, m, t, xi, nx, ny, g_pot=1, g_xi=1, filling_factor=1, epsilon=0.05):
        self.m = m
        self.t = t 
        self.nx = nx
        self.ny = ny
        self.g_pot = g_pot
        self.xi = xi
        self.ext_potential = np.zeros(self.nx)
        self.n_particles = 2*int(round(filling_factor*nx*ny))

        self.densities = {'up'  : np.zeros(2*self.nx).reshape(self.nx, 2),
                          'down': np.zeros(2*self.nx).reshape(self.nx, 2)}
        A,B = 0,1
        for spin in ('up', 'down'):
            self.densities[spin][:,B] = (np.ones(nx) + 
                                 epsilon*(np.random.rand(self.nx) - 0.5*np.ones(self.nx))) 

        self.figure_counter = 0


    def total_density(self):
        return np.sum(self.densities['up'], axis=1) + np.sum(self.densities['down'], axis=1)
    


    def hamiltonian(self, py, spin):
        H = np.zeros((2*self.nx)**2, dtype=complex).reshape(self.nx,2,self.nx,2)
        s = 1 if spin == 'up' else -1
        
        A,B = 0,1
        total_density = self.total_density()
        for i in xrange(self.nx):
           H[i,:,i,:] = np.array(
                    [self.xi + 1./self.m*(2 - math.cos(py)) 
                     + self.ext_potential[i]  
                     + self.g_pot*(total_density[i] - self.densities[spin][i, A]),
                    -2j*s*self.t*math.sin(py),
                    2j*s*self.t*math.sin(py), 
                    -self.xi - 1./self.m*(2 - math.cos(py)) 
                    + self.ext_potential[i]
                    + self.g_pot*(total_density[i] - self.densities[spin][i, B])]
                    ).reshape(2,2)
           if i > 0:
                H[i-1,:,i,:] = np.array([0.5/self.m, -1j*self.t,
                                         -1j*self.t, -0.5/self.m]).reshape(2,2)
                H[i,:,i-1,:] = np.array([0.5/self.m, 1j*self.t,
                                         1j*self.t, -0.5/self.m]).reshape(2,2)
        return H.reshape(2*self.nx,2*self.nx)


    def get_wf_density(self, wf, band_index):
        return np.sum(np.absolute(wf)**2, axis=band_index)
    
    
    def get_densities_and_efermi(self):
        pyrange = np.linspace(0, 2*math.pi, self.ny, endpoint=False)
    
        State = namedtuple('State', ['energy', 'wf', 'py', 'spin']) 
        states = []
        for py in pyrange:
            energies_up, wfs_up = eigh(self.hamiltonian(py, 'up'))
            energies_down, wfs_down = eigh(self.hamiltonian(py, 'down'))
    
            states_up = [State(*s) for s in izip(energies_up, wfs_up.T, 
                                                repeat(py), repeat('up'))]
            states_down = [State(*s) for s in izip(energies_down, wfs_down.T, 
                                                repeat(py), repeat('down'))]
    
            states.extend(states_up)
            states.extend(states_down)
    
        states.sort(key=lambda state: state.energy)

        densities = {'up'  : np.zeros(2*self.nx).reshape(self.nx, 2),
                     'down': np.zeros(2*self.nx).reshape(self.nx, 2)}
    
        for s in states[:self.n_particles]:
            densities[s.spin] += np.absolute(s.wf.reshape(self.nx,2))**2/self.ny
    
        fermi_energy = states[self.n_particles-1].energy
    
        return densities, fermi_energy


    def update_densities(self, sm_f):
        densities_new, fermi_energy = self.get_densities_and_efermi()
        for spin in ('up', 'down'): 
            self.densities[spin] = (sm_f*self.densities[spin] 
                                    + (1 - sm_f)*densities_new[spin])
        self.fermi_energy = fermi_energy
#        self.density['up'], self.density['down'], self.fermi_energy = \
#              (sm_f*self.density['up'] + (1 - sm_f)*d_up_new,
#               sm_f*self.density['down'] + (1 - sm_f)*d_down_new, 
#               fermi_energy)


    def hartree_fock_iterate(self, sm_f, n_iter):
        for i in xrange(n_iter):
            print '{}th Hartree-Fock iteration'.format(i)
            self.update_densities(sm_f)


#    def plot_densities(self):
#        density_figure = plt.figure(self.figure_counter)
#        self.figure_counter += 1
#
#        plt.plot(range(0,ham.nx), self.densities['up'],   color='blue')
#        plt.plot(range(0,ham.nx), self.densities['down'], color='green')
#    
#        plt.xlim(0,ham.nx-1)
#        density_figure.show()


    def plot_bands(self, spin):
        band_figure = plt.figure(self.figure_counter)
        self.figure_counter += 1
    
        energies = []
        Npy = 41
        pyrange = np.linspace(-math.pi, math.pi, Npy)
        for py in pyrange:
            energies.append(
                eigvalsh(self.hamiltonian(py, spin)))
        energies = np.array(energies).T
        
        for band in energies:
            plt.plot(pyrange, band, color='black')
    
        plt.plot(pyrange, self.fermi_energy*np.ones(len(pyrange)), color='red')
        plt.xlim(-math.pi, math.pi)
        
        band_figure.show()


if __name__ == '__main__':
    ham = TwoLevelStripeHam(m=1, t=0.4, xi=-1.5, 
                            nx=40, ny=60, g_pot=4, g_xi=0, filling_factor=0.98)
    ham.ext_potential = get_external_potential(ham.nx, length=6, strength=1.2)
#    ham.ext_xi = get_external_xi(ham.nx, 0.01, 1)

    ham.hartree_fock_iterate(sm_f=0., n_iter=5)
#    ham.plot_densities()
    ham.plot_bands('up')
    ham.plot_bands('down')
    raw_input()
