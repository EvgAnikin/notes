import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import bisect, fsolve


class Semiconductor:
    def __init__(self, gamma1, gamma2, gamma3, Ec, Ev):
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3
        self.Ec = Ec
        self.Ev = Ev
        self.beta = gamma1 + 2*gamma2


hg = Semiconductor(4.1, 0.5, 1.3, -0.303, 0)
cd = Semiconductor(1.47, -0.28, 0.03, 1.03, -0.57)
hbar_2m = 0.0381
Ep = 18.8


def get_eta(sc, E):
    eta_1 = 0.5/sc.beta*(
                -(2./3*Ep + sc.beta*(sc.Ec - E) - sc.Ev + E) + 
                math.sqrt((2./3*Ep + sc.beta*(sc.Ec - E) - sc.Ev + E)**2 + 
                          4*(sc.beta)*(sc.Ec - E)*(sc.Ev - E)))
    eta_2 = 0.5/sc.beta*(
                -(2./3*Ep + sc.beta*(sc.Ec - E) - sc.Ev + E) - 
                math.sqrt((2./3*Ep + sc.beta*(sc.Ec - E) - sc.Ev + E)**2 + 
                          4*(sc.beta)*(sc.Ec - E)*(sc.Ev - E)))
    return eta_1, eta_2


def normalize(row):
    return row/(math.sqrt(np.sum(row[:2]**2)))


def cd_row(eta, E, d):
    kappa = math.sqrt(-eta/hbar_2m)
    row = normalize(np.array([kappa, kappa**2, cd.Ec - E + eta, 
                    cd.beta*kappa*(cd.Ec - E + eta)]))
    return row


def hg_row(eta, E, d):
    if eta < 0:
        kappa = math.sqrt(-eta/hbar_2m)
        row = normalize(np.array([kappa*math.cosh(kappa*d/2),
                        -kappa**2*math.sinh(kappa*d/2),
                        -(hg.Ec - E - eta)*math.sinh(kappa*d/2),
                        hg.beta*kappa*math.cosh(kappa*d/2)*(hg.Ec - E + eta)]))
        return row
    elif eta > 0:
        k = math.sqrt(eta/hbar_2m)
        row = normalize(np.array([k*math.cos(k*d/2),
                         k**2*math.sin(k*d/2),
                        -(hg.Ec - E - eta)*math.sin(k*d/2),
                        hg.beta*k*math.cos(k*d/2)*(hg.Ec - E + eta)]))
        return row
    elif eta == 0:
        return hg_row(1e-7, E, d)


def determinant(E, d):
    eta_1_cd, eta_2_cd = get_eta(cd, E)
    eta_1_hg, eta_2_hg = get_eta(hg, E)

    mat = np.array([cd_row(eta_1_cd, E, d), cd_row(eta_2_cd, E, d), 
                    hg_row(eta_1_hg, E, d), hg_row(eta_2_hg, E, d)])
    return np.linalg.det(mat)
    

def plot_det(width):
    energies = np.linspace(cd.Ev, cd.Ec, 600)
    dets = [determinant(e, d=width) for e in energies]
    plt.plot(energies, dets)
    plt.show()


def plot_light_hole_level():
    widths = np.linspace(4,8,51)
    energies = []
    for width in widths:
        energies.append(bisect(lambda e: determinant(e, width), 
                               -0.2, 0.2))

    plt.plot(widths, energies)


def heavy_hole_levels(d):
    k_Ev = math.sqrt(-cd.Ev/hbar_2m)
    k_max = k_Ev/math.sqrt(hg.beta)
    n_roots = int(k_max*d/2./math.pi) + 1
    
    guesses = (np.linspace(0, n_roots-1, n_roots)*math.pi + math.pi/4)*2/d
    if guesses[-1] > k_max:
        epsilon = 1e-2
        guesses[-1] = k_max - epsilon
    
    equation = lambda k: k*np.tan(k*d/2) - (np.sqrt(
                           k_Ev**2*cd.beta/hg.beta**2 - 
                           cd.beta/hg.beta*k**2
                           ))
    
    levels= []
    for guess in guesses:
        k = fsolve(equation, guess)[0]
        levels.append(-hbar_2m*k**2*hg.beta)
    return np.array(levels)


def plot_heavy_hole_levels():
    for d in np.linspace(4,8,100):
        levels = heavy_hole_levels(d)
        plt.scatter([d]*len(levels), levels, marker='.')
    plt.show()


if __name__ == '__main__':
#    plot_det(5)
#    root = bisect(lambda e: determinant(e, d=4), -0.2,0.2)
#    print determinant(0, d=4)
#    print [math.log(abs(determinant(10**(-i), d=5))) for i in xrange(5)]
    plot_light_hole_level()
    plot_heavy_hole_levels()
    plt.show()
#    plot_det(5)
