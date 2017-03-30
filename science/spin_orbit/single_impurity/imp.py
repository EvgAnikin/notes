import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scint
from scipy.optimize import brentq

def hamiltonian(px, py, xi, m, t):
    return np.array([[xi + 1./m*(2 - math.cos(px) - math.cos(py)), 
                      2*t*(math.sin(px) - 1j*math.sin(py))],
                     [2*t*(math.sin(px) + 1j*math.sin(py)),
                      -xi - 1./m*(2 - math.cos(px) - math.cos(py))]])

def energy_square(px, py, xi, m, t):
    return ((xi + 1./m*(2 - math.cos(px) - math.cos(py)))**2
            + 4*t**2*(math.sin(px)**2 + math.sin(py)**2))
                    

def green_function(omega, px, py, xi, m, t):
    return ((np.identity(2)*omega + hamiltonian(px,py,xi,m,t))
            /(omega**2 - energy_square(px,py,xi,m,t)))


def green_function_11(omega, xi, m, t):
    g00_11 = lambda px, py:\
             ((omega + xi + 1./m*(2 - math.cos(px) - math.cos(py)))
              /(omega**2 - energy_square(px, py, xi, m, t))/(2*math.pi)**2)
    return scint.dblquad(g00_11, -math.pi, math.pi, 
                                 lambda x: -math.pi, lambda x: math.pi)[0]


def gf_divergent_p(omega, xi, m, t, p_max):
    g00_11 = lambda px, py:\
             ((omega + xi + 1./m*(2 - math.cos(px) - math.cos(py)))
              /(omega**2 - energy_square(px, py, xi, m, t))/(2*math.pi)**2)
    return scint.dblquad(g00_11, -p_max, p_max, 
                         lambda x: -math.sqrt(p_max**2 - x**2), 
                         lambda x: math.sqrt(p_max**2 - x**2))[0]


def gf_regular_p(omega, xi, m, t, p_max):
    g00_11 = lambda px, py:\
             ((omega + xi + 1./m*(2 - math.cos(px) - math.cos(py)))
              /(omega**2 - energy_square(px, py, xi, m, t)))/(2*math.pi)**2
    return 2*scint.dblquad(g00_11, -math.pi, math.pi, 
                         lambda x: 0 if abs(x) > p_max else math.sqrt(p_max**2 - x**2), 
                         lambda x: math.pi)[0]


def gf_divergent_p_approx(omega, xi, m, t, p_max):
    return -(1./(8*math.pi)/(m*(4*t**2 + xi/m))*
            (p_max**2 + (2*m*(omega + xi) - (xi**2 - omega**2)/(4*t**2 + xi/m))*
                        math.log(1 + (4*t**2 + xi/m)*p_max**2/(xi**2 - omega**2))))


def gf_full(omega, xi, m, t, p_max):
    return (gf_divergent_p_approx(omega, xi, m, t, p_max) +
            gf_regular_p(omega, xi, m, t, p_max))


def det_green_function(omega, xi, m, t, dE):
    g00_11 = lambda px, py:\
             ((omega + xi + 1./m*(2 - math.cos(px) - math.cos(py)))
              /(omega**2 - energy_square(px, py, xi, m, t)))
    g00_22 = lambda px, py:\
             ((omega - xi - 1./m*(2 - math.cos(px) - math.cos(py)))
              /(omega**2 - energy_square(px, py, xi, m, t)))
    return ((1 - dE*scint.dblquad(g00_00, -math.pi, math.pi, 
                                 lambda x: -math.pi, lambda x: math.pi)[0])
            *(1 - dE*scint.dblquad(g00_11, -math.pi, math.pi, 
                                 lambda x: -math.pi, lambda x: math.pi)[0]))
 
                
if __name__ == '__main__':
    epsilon = 1e-5
    xi, m, t = -0.03, 0.1, 0.5
    p_max = 0.05
    max_omega = (1 - epsilon)*math.sqrt(energy_square(0, 0, xi, m, t))*(1 - epsilon)

    N = 10
#    omegas = np.linspace(-max_omega, max_omega, N)
    omegas = 2*max_omega*np.power(10, np.linspace(-5,0,N)) - max_omega
#    gf = [gf_divergent_p_approx(xi + domega, xi, m, t, p_max) + gf_regular 
#             for domega in domegas]

    gf = []
    counter = 0
    for omega in omegas:
        counter += 1
        print '{}th value'.format(counter)
        gf.append(green_function_11(omega, xi, m, t))
    
    gf = np.array(gf)

    output = open('gf_array', 'w')
    np.savez(output, omegas=omegas, gf=gf)

    plt.plot(omegas, gf)
    plt.plot(omegas[::-1], -gf[::-1])
    plt.show()

#    dets = [det_green_function(omega, xi, m, t, dE) for omega in omegas]
#    print brentq(det_green_function, 0, max_omega, args = (xi, m, t, dE))

#    print det_green_function(0, xi, m, t, dE)/dE**2
#    print det_green_function(max_omega, xi, m, t, dE)/dE**2
