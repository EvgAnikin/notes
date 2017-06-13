import kwant
import matplotlib.pyplot as plt
import math
import numpy as np


def qpc_shape(position, a, b, sigma, x_max):
    x,y = position
    return abs(y) < a - (a - b)*math.exp(-x**2/2./sigma**2) and abs(x) < x_max


def build_lead(a, x_max, direction, sqlat):
    lead = kwant.Builder(kwant.TranslationalSymmetry([direction, 0]))
    lead[(sqlat(direction*(x_max-1), y) for y in range(-a+1, a))] = 4
    lead[sqlat.neighbors()] = -1
    return lead


def build_bar(nx, ny, xi, m, t):
    bar = kwant.Builder()
    sqlat = kwant.lattice.square()
    bar_shape = lambda x,y: abs(x) < nx and abs(y) < ny
    bar[sqlat.shape(bar_shape, (0,0))] = tinyarray.array([[xi + 2/m, 0],
                                         [0, -xi - 2/m]])

    bar[sqlat.neighbors()] = -1

    return qpc


if __name__ == '__main__':
    qpc = build_qpc(a=15, b=4, sigma=6, x_max=20) 
    kwant.plot(qpc)

    energies = np.linspace(0,1,100)
    conductances = [kwant.smatrix(qpc, en).transmission(1, 0)
                                            for en in energies]
    plt.plot(energies, conductances)
    plt.show()
