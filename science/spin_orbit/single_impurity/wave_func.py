import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scint
from scipy.optimize import brentq


def energy_square(px, py, xi, m, t):
    return ((xi + 1./m*(2 - math.cos(px) - math.cos(py)))**2
            + 4*t**2*(math.sin(px)**2 + math.sin(py)**2))


def green_function_11(omega, xi, m, t, n):
    g00_11 = lambda px, py:\
             ((omega + xi + 1./m*(2 - math.cos(px) - math.cos(py)))*math.cos(px*n)
              /(omega**2 - energy_square(px, py, xi, m, t))/(2*math.pi)**2)
    return scint.dblquad(g00_11, -math.pi, math.pi, 
                                 lambda x: -math.pi, lambda x: math.pi)[0]

xi,m,t = -0.1, 1, 0.5
omega = -0.09

gf = []
for n in xrange(0,10):
   gf.append(green_function_11(omega, xi, m, t, n))

plt.plot(gf)
plt.show()
