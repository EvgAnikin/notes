import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

hbar = 1.0545718e-27
me = 9.10938356e-28
eV = 1.6021766209e-12
a = 6.46e-8
nm = 1e-7

gamma_1_hgte = 4.1
gamma_2_hgte = -0.28
gamma_3_hgte = 1.3

gamma_1_cdte = 1.47
gamma_2_cdte = -0.28
gamma_3_cdte = 0.03

Ev = 0.57

def heavy_hole_levels(d):
    k_Ev = math.sqrt(2*me*Ev*eV/hbar**2*nm**2)
    k_max = k_Ev/math.sqrt((gamma_1_hgte - 2*gamma_2_hgte))
    n_roots = int(k_max*d/2./math.pi) + 1
    
    guesses = (np.linspace(0, n_roots-1, n_roots)*math.pi + math.pi/4)*2/d
    if guesses[-1] > k_max:
        epsilon = 1e-2
        guesses[-1] = k_max - epsilon
    
    equation = lambda k: k*np.tan(k*d/2) - (np.sqrt(
            k_Ev**2*(gamma_1_cdte - 2*gamma_2_cdte)/(gamma_1_hgte - 2*gamma_2_hgte)**2 - 
            (gamma_1_cdte - 2*gamma_2_cdte)/(gamma_1_hgte - 2*gamma_2_hgte)*k**2
            ))
#    ks = np.linspace(0, k_max, 100)
#    plt.plot(ks/math.pi, equation(ks))
#    plt.scatter(guesses/math.pi, np.zeros(len(guesses))) 
#    plt.ylim(-5,5)
#    plt.show()
    
    levels= []
    for guess in guesses:
        k = fsolve(equation, guess)[0]
        levels.append(-hbar**2/2/me/nm**2/eV*k**2*(gamma_1_hgte - 2*gamma_2_hgte))
    return np.array(levels)
    

if __name__ == '__main__':
    for d in np.linspace(4,8,100):
        levels = heavy_hole_levels(d)
        plt.scatter([d]*len(levels), levels, marker='.')
    plt.show()
