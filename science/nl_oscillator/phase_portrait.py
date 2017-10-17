#encoding:utf-8
from __future__ import print_function
from __future__ import division
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def osc_energy(re, im, delta, beta, E):
    return delta*(re**2 + im**2) + beta/2*(re**2 + im**2)**2 + 2*E*re


N = 200

re = np.linspace(-30, 20, N)
im = np.linspace(-20, 20, N)


delta = -1
beta = 1/12**2
E = 2/3/math.sqrt(3)*abs(delta)**1.5/beta**0.5 * 0.8
print(E)

root_sep = fsolve(lambda x: x**3 + delta/beta*x + E/beta, [math.sqrt(abs(delta)/beta)])[0]
separatrix_energy = osc_energy(root_sep, 0, delta, beta, E)

root_min = fsolve(lambda x: x**3 + delta/beta*x + E/beta, [-math.sqrt(abs(delta)/beta)])[0]
minimum_energy = osc_energy(root_min, 0, delta, beta, E)

re, im = np.meshgrid(re, im)
energy = osc_energy(re, im, delta, beta, E)

matplotlib.rcParams.update({'font.size': 25})

n_lines_below = 8
plt.contour(re, im, energy, np.linspace(minimum_energy, separatrix_energy, n_lines_below))

delta_energy = (separatrix_energy - minimum_energy)/(n_lines_below - 1)
n_lines_above = 190
plt.contour(re, im, energy, [separatrix_energy + k*delta_energy for k in xrange(2,n_lines_above, 5)])

plt.show()
plt.savefig('figure.png')
