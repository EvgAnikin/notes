import kwant
import matplotlib.pyplot as plt
import math
import numpy as np


def build_lead(a, x_max, direction, sqlat):
    lead = kwant.Builder(kwant.TranslationalSymmetry([direction, 0]))
    lead[(sqlat(direction*(x_max-1), y) for y in range(-a+1, a))] = 4
    lead[sqlat.neighbors()] = -1
    return lead


def build_disordered_wire(a, disorder, x_max):
    wire = kwant.Builder()
    sqlat = kwant.lattice.square()

    for x in range(-x_max+1, x_max):
        for y in range(-a+1, a):
            sigma = 0.5*x_max
            wire[sqlat(x,y)] = 4 + disorder*np.random.rand()*math.exp(-0.5*x**2/sigma**2)

    wire[sqlat.neighbors()] = -1

    left_lead  = build_lead(a, x_max, -1, sqlat)
    right_lead = build_lead(a, x_max,  1, sqlat)
    wire.attach_lead(left_lead)
    wire.attach_lead(right_lead)
    wire = wire.finalized()
    return wire


if __name__ == '__main__':
    energy = 0.3
    length = list(range(1,1001,100))
    conductances = []
    for l in length:
        wire = build_disordered_wire(a=40, disorder=0.3, x_max=l) 
        conductances.append(kwant.smatrix(wire, energy).transmission(1, 0))

    resistances = 1/np.array(conductances)
    plt.plot(length, resistances)
    plt.ylim(0, 1.1*np.amax(resistances))
    plt.show()
