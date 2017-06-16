import math
import matplotlib.pyplot as plt
import numpy as np


def plot_histogram(energies, bins):
    values, base = np.histogram(energies, bins=bins)
    emax, emin = np.amax(energies), np.amin(energies)
    density = values/float(energies.size)*bins/(emax - emin)
    plt.plot((base[:-1] + base[1:])/2, density)
    plt.ylim(0, np.amax(density)*1.2)
    plt.show()


def show_top_ins_dos(xi, m, t):
    def top_ins_energy(px,py):
        return np.sqrt( (xi + 1./m*(2 - np.cos(px) - np.cos(py)))**2 + 
                         4*t**2*(np.sin(px)**2 + np.sin(py)**2) )
    NX, NY = 1000, 1000
    px, py = np.meshgrid(np.linspace(0, 2*math.pi, NX, endpoint=False),
                         np.linspace(0, 2*math.pi, NY, endpoint=False))
    energies = top_ins_energy(px,py)
    plot_histogram(energies, 100)


if __name__ == '__main__':
    xi, m, t = 0.1, 1, 0.3
    show_top_ins_dos(xi, m, t)
