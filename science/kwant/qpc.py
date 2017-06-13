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


def build_qpc(a, b, sigma, x_max):
    qpc = kwant.Builder()
    sqlat = kwant.lattice.square()
    qpc[sqlat.shape(lambda pos: qpc_shape(pos, a, b, sigma, x_max), (0,0))] = 4
    qpc[sqlat.neighbors()] = -1

    left_lead  = build_lead(a, x_max, -1, sqlat)
    right_lead = build_lead(a, x_max,  1, sqlat)
    qpc.attach_lead(left_lead)
    qpc.attach_lead(right_lead)
    qpc = qpc.finalized()
    return qpc


if __name__ == '__main__':
    qpc = build_qpc(a=15, b=4, sigma=6, x_max=20) 
    kwant.plot(qpc)

    energies = np.linspace(0,1,100)
    conductances = [kwant.smatrix(qpc, en).transmission(1, 0)
                                            for en in energies]
    plt.plot(energies, conductances)
    plt.show()
