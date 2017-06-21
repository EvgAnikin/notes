from functools import partial
import kwant
import matplotlib.pyplot as plt
import math
import numpy as np
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



def build_bar(nx, ny, xi, m, t, E_lead, t_lead, cylinder=False):
    bar = kwant.Builder()

    lat = kwant.lattice.general([(1,0), (0,1)], [(0,0), (0,0)])
    a,b = lat.sublattices 
    bar_shape = lambda pos: abs(pos[0]) < nx and abs(pos[1]) < ny

    def onsite_with_disorder(sublattice, site, disorder_strength=0, seed=''):
        mult = 1 if sublattice == 'a' else -1
        return  (mult*(xi + 2/m) + 
                      disorder_strength*(2*kwant.digest.uniform(repr(site.pos) + seed) - 1))

    U = 3
    disorder_strength = 0.7
    def onsite_with_impurities(sublattice, site, probability=0, seed=''):
        rand_num = kwant.digest.uniform(repr(site.pos) + seed)
        mult = 1 if sublattice == 'a' else -1
#        if rand_num < probability:
#            return mult*(-xi + 2./m)
#        else:
#            return mult*(xi + 2./m)
        if rand_num < probability:
            return  (mult*(xi + 2./m) + U + 
                     disorder_strength*(2*kwant.digest.uniform(repr(site.pos) + seed) - 1))
        else:
            return mult*(xi + 2./m)


    bar[a.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'a')
    bar[b.shape(bar_shape, (0,0))] = partial(onsite_with_impurities, 'b')

#    bar[a.shape(bar_shape, (0,0))] = xi + 2./m
#    bar[b.shape(bar_shape, (0,0))] = -xi - 2./m

#    for i in range(4, nx):
#        bar[a(0, i)] = 100
#        bar[b(0, i)] = 100
#
#        bar[a(0,-i)] = 100
#        bar[b(0,-i)] = 100

#    bar[a(0,0)] += U
#    bar[b(0,0)] += U

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
    bar.attach_lead(left_lead, add_cells=1)
    bar[((a(-nx, i), b(-nx+1,i)) for i in range(-ny+1,ny))] = -t_lead

    right_lead = build_lead(ny, t_lead, E_lead, 1, a)
    bar.attach_lead(right_lead, add_cells=1)
    bar[((a(nx, i), b(nx-1,i)) for i in range(-ny+1,ny))] = -t_lead

    bar = bar.finalized()
    return bar


if __name__ == '__main__':
    xi, m, t = -0.1, 1, 0.4
    E_lead, m_lead = -0.8, 0.5
    bar = build_bar(50, 30, xi, m, t, E_lead, m_lead, cylinder=False)

    energy = -0.04
    probability = 0.005
    seed = repr(np.random.rand())
    smatrix = kwant.smatrix(bar, energy, args = [probability, seed])
    print('seed : {}, transmission: {}'.format(seed, smatrix.transmission(1, 0)))

    local_dos = kwant.ldos(bar, energy=energy, args=[probability, seed])
    kwant.plot(bar, #site_color = -local_dos,
               site_color=np.minimum(-np.log(local_dos), 8*np.ones_like(local_dos)),
               hop_color='white', cmap='afmhot')


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

#    emin, emax = -2*abs(xi), 2*abs(xi)
#    probability = 0.015
#
#    energies = np.linspace(emin, emax, 16)
#    conductances = []
#    counter = 0
#    for en in energies:
#        print(counter)
#        counter += 1
#        conductances.append(kwant.smatrix(bar, en, args=[probability]).transmission(1, 0))
#
#    plt.plot(energies, conductances)
#    plt.xlim(emin, emax)
#    plt.ylim(0,2)
#    plt.show()
