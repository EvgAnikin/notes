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
    def __init__(self, U, nx, ny): 
        self.U = U
        self.nx = nx
        self.ny = ny

    
    def get_impurity_potential(self, pos, disorder, probability, seed=''):
        rand_num_1 = kwant.digest.uniform(repr(pos) + seed + '1')
        rand_num_2 = kwant.digest.uniform(repr(pos) + seed + '2')
        x,y = pos
    
        if rand_num_1 > probability or abs(x) >= self.nx or abs(y) >= self.ny:
            return 0
        else:
            return (self.U + disorder*(2*rand_num_2 - 1))

    
    def get_impurities_array(self, disorder, probability, seed=''):
        energies = []
        for x in range(-self.nx+1, self.nx):
            for y in range(-self.ny+1, self.ny):
                e = self.get_impurity_potential((x, y), disorder, probability, seed)
                if e != 0:
                    energies.append(e)
        return energies


def build_bar(nx, ny, xi, m1, m2, t, E_lead, t_lead, 
              dis_gen=None, cylinder=False, with_leads=True):
    bar = kwant.Builder()

    lat = kwant.lattice.general([(1,0), (0,1)], [(0,0), (0,0)])
    a,b = lat.sublattices 
    bar_shape = lambda pos: abs(pos[0]) < nx and abs(pos[1]) < ny

    def onsite_with_impurities(sublattice, site, disorder=0, probability=0, seed=''):
        pos = (int(site.pos[0]), int(site.pos[1]))
        dis = dis_gen.get_impurity_potential(pos, disorder, probability, seed)
        if sublattice == 'a':
            return ( xi + 2./m1) + dis
        else:
            return (-xi - 2./m2) + dis

    if dis_gen:
        bar[a.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'a')
        bar[b.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'b')
    else:
        bar[a.shape(bar_shape, (0,0))] = ( xi + 2./m1)
        bar[b.shape(bar_shape, (0,0))] = (-xi - 2./m2)
        

    bar[kwant.builder.HoppingKind((-1,0), a, a)] = -0.5/m1
    bar[kwant.builder.HoppingKind((-1,0), a, b)] = -1j*t
    bar[kwant.builder.HoppingKind((-1,0), b, a)] = -1j*t
    bar[kwant.builder.HoppingKind((-1,0), b, b)] = 0.5/m2

    bar[kwant.builder.HoppingKind((0,-1), a, a)] = -0.5/m1
    bar[kwant.builder.HoppingKind((0,-1), a, b)] = -t
    bar[kwant.builder.HoppingKind((0,-1), b, a)] = t
    bar[kwant.builder.HoppingKind((0,-1), b, b)] = 0.5/m2

    if cylinder:
        for i in range(-nx+1, nx):
            bar[a(i, ny-1), a(i,-ny+1)] = -0.5/m1
            bar[a(i, ny-1), b(i,-ny+1)] = -t
            bar[b(i, ny-1), a(i,-ny+1)] = t
            bar[b(i, ny-1), b(i,-ny+1)] = 0.5/m2

    if with_leads:
        left_lead = build_lead(ny, t_lead, E_lead, -1, a, cylinder)
        bar.attach_lead(left_lead, add_cells=1)
        bar[((a(-nx, i), b(-nx+1,i)) for i in range(-ny+1,ny))] = -t_lead
    
        right_lead = build_lead(ny, t_lead, E_lead, 1, a)
        bar.attach_lead(right_lead, add_cells=1)
        bar[((a(nx, i), b(nx-1,i)) for i in range(-ny+1,ny))] = -t_lead
    else:
        pass
        for i in range(-nx+1, nx):
            bar[a(i, ny-1), a(i,-ny+1)] = -0.5/m1
            bar[a(i, ny-1), b(i,-ny+1)] = -t
            bar[b(i, ny-1), a(i,-ny+1)] = t
            bar[b(i, ny-1), b(i,-ny+1)] = 0.5/m2

        for i in range(-ny+1, ny):
            bar[a(nx-1, i), a(-nx+1, i)] = -0.5/m1
            bar[a(nx-1, i), b(-nx+1, i)] = -1j*t 
            bar[b(nx-1, i), a(-nx+1, i)] = -1j*t 
            bar[b(nx-1, i), b(-nx+1, i)] = 0.5/m2

    bar = bar.finalized()
    return bar


def get_level_of_impurity_function():
    input_file = open('gf_array', 'rb')
    np_file = np.load(input_file)
    
    omegas = np_file['omegas']
    gf = np_file['gf']
    return interp1d(1/gf, omegas, kind='cubic'), (np.amin(1/gf), np.amax(1/gf))


def plot_levels(levels, color, ymin=0, ymax=0.5):
    for e in levels:
        plt.plot([e,e], [ymin, ymax], color=color)


def plot_conductance_vs_energy(bar, disorder, probability, seed):
    emin, emax = -50, 50
    energies = np.linspace(emin, emax, 30)
    conductances = []
    counter = 0
    for en in energies:
        print(counter)
        counter += 1
        conductances.append(kwant.smatrix(bar, en, 
                            args=[disorder, probability, seed]).transmission(1, 0))
                                                         
    plt.plot(energies, conductances)
    plt.xlim(emin, emax)
    plt.show()


def plot_conductance_vs_disorder(bar, energy):
    probabilities = np.linspace(0, 0.1, 20)
    disorder = 500
    conductances = []
    errors = []
    
    for counter, prob in enumerate(probabilities):
        print(counter)
        cond = []
        n_sim = 100
        for i in range(n_sim):
            seed = str(np.random.rand())
            smatrix = kwant.smatrix(bar, energy, args = [disorder, prob, seed])
            cond.append(smatrix.transmission(1, 0))
        conductances.append(np.average(cond))
        errors.append(np.std(cond))

    ax = plt.subplot(111)
#    ax.set_xscale('log')
    ax.errorbar(probabilities, conductances, yerr=errors)
#    ax.set_xlim(xmin=1)
    ax.set_ylim(ymin=-0.1, ymax=1.1)
    plt.show()

    out = open('conductances_plot', 'wb')
    np.savez(out, probabilities=probabilities, conductances=conductances, errors=errors)


def plot_local_dos(bar, energy, disorder, probability, seed):
    smatrix = kwant.smatrix(bar, energy, args = [disorder, probability, seed])
    print('seed : {}, transmission: {}'.format(seed, smatrix.transmission(1, 0)))

    local_dos = kwant.ldos(bar, energy=energy, args=[disorder, probability, seed])
    kwant.plot(bar, #site_color = -local_dos,
               site_color=np.minimum(-np.log(local_dos), 10*np.ones_like(local_dos)),
               hop_color='white', cmap='afmhot')


if __name__ == '__main__':
    nx, ny = 50, 100 
    a = 5
    A, B, D, M = 364.5, -686, 0, -10
    
    xi = M
    t = A/a
    m1 = -0.25*a**2/(B + D)
    m2 = -0.25*a**2/(B - D)

    E_lead, m_lead = -150, 0.2*m1
#    xi, m1, m2, t = -0.02, 1, 1, 0.4

    dis_gen = DisorderGenerator(U=50, nx=nx, ny=ny)
    bar = build_bar(nx, nx, xi, m1, m2, t, E_lead, m_lead, 
                   dis_gen=dis_gen)
    bar_without_leads = build_bar(nx, nx, xi, m1, m2, t, E_lead, m_lead, 
                   dis_gen=dis_gen, with_leads=False)

    probability = 0.05
    disorder = 200

#    plot_local_dos(bar, 0, 0, 0, seed)
    plot_conductance_vs_energy(bar, disorder, probability, seed='')
#    plot_conductance_vs_disorder(bar, energy=0)

#    imp_energies = np.array(dis_gen.get_impurities_array(probability, seed))
#    get_imp_level, (low, high) = get_level_of_impurity_function()
#    imp_levels = get_imp_level(np.array([e for e in imp_energies if e > low and e < high]))

#    ham = bar_without_leads.hamiltonian_submatrix(args=[probability, seed])
#    energies = eigvalsh(ham)

#    plot_levels(imp_levels, color='green')
#    plot_levels(energies, color='blue')

