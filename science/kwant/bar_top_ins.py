from functools import partial
import kwant
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import eigvalsh
import sys
import tinyarray


def build_lead(ny, m, E, direction, lat, cylinder=False):
    lead = kwant.Builder(kwant.TranslationalSymmetry([direction, 0]))
    lead_shape = lambda pos: abs(pos[0]) < nx
    lead[(lat(0,i) for i in range(-ny+1,ny))] = 2/m + E
    lead[lat.neighbors()] = -0.5/m
    if cylinder:
        lead[lat(0,ny-1), lat(0,-ny+1)] = -0.5/m
    return lead


class DisorderGenerator:
    def __init__(self, U, disorder, nx, ny): 
        self.U = U
        self.disorder = disorder
        self.nx = nx
        self.ny = ny

    
    def get_impurity_potential(self, pos, probability, seed=''):
        rand_num_1 = kwant.digest.uniform(repr(pos) + seed + '1')
        rand_num_2 = kwant.digest.uniform(repr(pos) + seed + '2')
        x,y = pos
    
        if rand_num_1 > probability or abs(x) >= self.nx or abs(y) >= self.ny:
            return 0
        else:
            return (self.U + self.disorder*(2*rand_num_2 - 1))

    
    def get_impurities_array(self, probability, seed=''):
        energies = []
        for x in range(-self.nx+1, self.nx):
            for y in range(-self.ny+1, self.ny):
                e = self.get_impurity_potential((x, y), probability, seed)
                if e != 0:
                    energies.append(e)
        return energies


def build_bar(nx, ny, xi, m, t, E_lead, t_lead, cylinder=False):
    bar = kwant.Builder()

    lat = kwant.lattice.general([(1,0), (0,1)], [(0,0), (0,0)])
    a,b = lat.sublattices 
    bar_shape = lambda pos: abs(pos[0]) < nx and abs(pos[1]) < ny

    dis_gen = DisorderGenerator(U=-2.5, disorder=0.02, nx=nx, ny=ny)
    def onsite_with_impurities(sublattice, site, probability=0, seed=''):
        mult = 1 if sublattice == 'a' else -1
        return mult*(xi + 2./m) + dis_gen.get_impurity_potential(site.pos, probability, seed)


    bar[a.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'a')
    bar[b.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'b')

    bar[kwant.builder.HoppingKind((-1,0), a, a)] = -0.5/m
    bar[kwant.builder.HoppingKind((-1,0), a, b)] = -1j*t
    bar[kwant.builder.HoppingKind((-1,0), b, a)] = -1j*t
    bar[kwant.builder.HoppingKind((-1,0), b, b)] = 0.5/m

    bar[kwant.builder.HoppingKind((0,-1), a, a)] = -0.5/m
    bar[kwant.builder.HoppingKind((0,-1), a, b)] = -t
    bar[kwant.builder.HoppingKind((0,-1), b, a)] = t
    bar[kwant.builder.HoppingKind((0,-1), b, b)] = 0.5/m

    if cylinder:
        for i in range(-nx+1, nx):
            bar[a(i, ny-1), a(i,-ny+1)] = -0.5/m
            bar[a(i, ny-1), b(i,-ny+1)] = -t
            bar[b(i, ny-1), a(i,-ny+1)] = t
            bar[b(i, ny-1), b(i,-ny+1)] = 0.5/m

    left_lead = build_lead(ny, t_lead, E_lead, -1, a, cylinder)
    bar.attach_lead(left_lead, add_cells=0)
#    bar[((a(-nx, i), b(-nx+1,i)) for i in range(-ny+1,ny))] = -t_lead

    right_lead = build_lead(ny, t_lead, E_lead, 1, a)
    bar.attach_lead(right_lead, add_cells=0)
#    bar[((a(nx, i), b(nx-1,i)) for i in range(-ny+1,ny))] = -t_lead

    bar = bar.finalized()
    return bar, dis_gen


def get_level_of_impurity_function():
    input_file = open('gf_array', 'rb')
    np_file = np.load(input_file)
    
    omegas = np_file['omegas']
    gf = np_file['gf']
    return interp1d(1/gf, omegas, kind='cubic'), (np.amin(1/gf), np.amax(1/gf))


def plot_levels(levels, ymin, ymax, color):
    for e in levels:
        plt.plot([e,e], [ymin, ymax], color=color)


if __name__ == '__main__':
    xi, m, t = -0.1, 1, 0.4
    E_lead, m_lead = -0.8, 0.5
    bar, dis_gen = build_bar(10, 10, xi, m, t, E_lead, 
                                                m_lead, cylinder=False)

#    emin, emax = -0.01, 0.01
    probability = 0#1./(20*15)
    seed = repr(np.random.rand())

    imp_energies = np.array(dis_gen.get_impurities_array(probability, seed))
    get_imp_level, (low, high) = get_level_of_impurity_function()
    imp_levels = get_imp_level(np.array([e for e in imp_energies if e > low and e < high]))

    bar_hamiltonian = bar.hamiltonian_submatrix(args=[probability, seed])
    system_energies_in_gap = [e for e in eigvalsh(bar_hamiltonian) if e > -1.5*abs(xi) and
                                                                   e <  1.5*abs(xi)]
#
#    energies = np.linspace(emin, emax, 20)
#    conductances = []
#    counter = 0
#    for en in energies:
#        print(counter)
#        counter += 1
#        conductances.append(kwant.smatrix(bar, en, args=[probability, 
#                                                         seed]).transmission(1, 0))

    plot_levels(imp_levels, 0., 0.5, color='green')
    plot_levels(system_energies_in_gap, 0., 0.5, color='blue')
#    plt.plot(energies, conductances)
#    plt.xlim(emin, emax)
    plt.ylim(0,1)
    plt.show()
    
#    energy = 0.0
#    smatrix = kwant.smatrix(bar, energy, args = [probability, seed])
#    print('seed : {}, transmission: {}'.format(seed, smatrix.transmission(1, 0)))
#
#    local_dos = kwant.ldos(bar, energy=energy, args=[probability, seed])
#    kwant.plot(bar, #site_color = -local_dos,
#               site_color=np.minimum(-np.log(local_dos), 8*np.ones_like(local_dos)),
#               hop_color='white', cmap='afmhot')
#

#    concentrations = np.linspace(0, 1, 40)
#    conductances = []
#    
#    for counter, c in enumerate(concentrations):
#        print(counter)
#        smatrix = kwant.smatrix(bar, 0., args = [c, repr(np.random.rand())])
#        conductances.append(smatrix.transmission(1, 0))
#
#    plt.plot(concentrations, conductances)
##    plt.ylim(-0.1, 1.1)
#    plt.show()

