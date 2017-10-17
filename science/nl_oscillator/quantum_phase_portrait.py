#encoding:utf-8
from __future__ import print_function
from __future__ import division
import cmath
import math
import mpmath
import numpy as np
from numpy.linalg import eig, eigh, eigvalsh, lstsq
import matplotlib
import matplotlib.pyplot as plt


def nonlinear_oscillator(omega, beta, E, N):
    ham = np.zeros(N**2).reshape(N,N)

    for i in xrange(N):
        ham[i,i] = omega*i + beta/2*(i**2 - 1)
        if i < N - 1:
            ham[i,i+1] = math.sqrt(i+1)*E
            ham[i+1,i] = math.sqrt(i+1)*E

    return ham


def coherent_representation(state, zbar):
    polynom = 0
    monom = 1
    for i, c in enumerate(state):
        polynom += c*monom
        monom *= zbar/math.sqrt(i+1)

    return math.exp(-0.5*abs(zbar)**2)*polynom


N_E = 400
delta = -1
beta = 1/12**2
E = 2/3/math.sqrt(3)*abs(delta)**1.5/beta**0.5 * 0.8
print(E)

energies, vectors = eigh(nonlinear_oscillator(delta, beta, E, 500))

def plot_state(k, rerange=[-30, 20], imrange=[-20, 20], points=10000):
    mpmath.cplot(lambda z: coherent_representation(vectors[:, k], z), 
                        rerange, imrange, points)

plot_state(10)
