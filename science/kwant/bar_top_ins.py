from functools import partial
import kwant
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import tinyarray


def build_lead(ny, m, E, direction, lat):
    lead = kwant.Builder(kwant.TranslationalSymmetry([direction, 0]))
    lead_shape = lambda pos: abs(pos[0]) < nx
    lead[(lat(0,i) for i in range(-ny+1,ny))] = 2/m + E
    lead[lat.neighbors()] = -0.5/m
    return lead



def build_bar(nx, ny, xi, m, t, E_lead, t_lead):
    bar = kwant.Builder()

    lat = kwant.lattice.general([(1,0), (0,1)], [(0,0), (0,0)])
    a,b = lat.sublattices 
    bar_shape = lambda pos: abs(pos[0]) < nx and abs(pos[1]) < ny

    def onsite_with_disorder(sublattice, site, disorder_strength=0, seed=''):
        mult = 1 if sublattice == 'a' else -1
        return  (mult*(xi + 2/m) + 
                      disorder_strength*(2*kwant.digest.uniform(repr(site.pos) + seed) - 1))

    U = 2.5
    disorder_strength = 0.1
    def onsite_with_impurities(sublattice, site, probability=0, seed=''):
        rand_num = kwant.digest.uniform(repr(site.pos) + seed)
        mult= 1 if sublattice == 'a' else -1
        if rand_num < probability:
            return  (mult*(xi + 2/m) + U + 
                     disorder_strength*(2*kwant.digest.uniform(repr(site.pos) + seed) - 1))
        else:
            return mult*(xi + 2/m)


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

    left_lead = build_lead(ny, t_lead, E_lead, -1, a)
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
    bar = build_bar(80, 40, xi, m, t, E_lead, m_lead)

#    smatrix = kwant.smatrix(bar, 0., args = [0.02])
#    print(smatrix.transmission(1, 0))
#
#    local_dos = kwant.ldos(bar, energy=0.00, args=[0.02])
#    kwant.plot(bar, site_color=np.minimum(-np.log(local_dos), 8*np.ones_like(local_dos)),
#               hop_color='white', cmap='afmhot')
#

#    concentrations = np.linspace(0, 0.1, 10)
#    conductances = []
#    
#    for counter, c in enumerate(concentrations):
#        print(counter)
#        smatrix = kwant.smatrix(bar, 0., args = [c, repr(np.random.rand())])
#        conductances.append(smatrix.transmission(1, 0))
#
#    plt.plot(concentrations, conductances)
#    plt.show()

    energies = np.linspace(-abs(xi), abs(xi), 60)
    conductances = [kwant.smatrix(bar, en).transmission(1, 0)
                                            for en in energies]
    plt.plot(energies, conductances)
    plt.show()
