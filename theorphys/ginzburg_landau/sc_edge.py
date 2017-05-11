import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

def edge_free_energy(fields, N, kappa, dx):
    # Setting constraints
    fields_reshaped = fields.reshape(2,2*N)
    psi, ay = fields_reshaped[0], fields_reshaped[1]

    psi = np.insert(np.append(psi, 1), 0, 0)
    ay = np.insert(np.append(ay, 0), 0, -N*dx/math.sqrt(2))

    part_without_ders = (1/kappa**2*np.sum(ay**2*psi**2)
                                  - np.sum(psi**2) + 
                                    np.sum(psi**4)/2)*dx
    part_with_ders = np.sum(np.diff(psi)**2)/dx + \
                    np.sum((np.diff(ay) - 1/math.sqrt(2)*dx*np.ones(2*N + 1))**2)/dx
    return part_without_ders + part_with_ders


def free_energy_density(psi, ay, N, kappa, dx):
    return ((1/kappa**2*ay**2*psi**2 - psi**2 + psi**4/2)[:-1] 
            + np.diff(psi)**2/dx**2 
            + (np.diff(ay)/dx - 1/math.sqrt(2)*np.ones(2*N - 1))**2)
    


def initial_guess(N, dx):
    ay = np.zeros(2*N)
    ay[0:N] = np.linspace(-N*dx, 0, N, endpoint=False)
    psi = np.zeros(2*N)
    psi[N:] = np.ones(N)
    return np.array([psi, ay]).reshape(4*N)


def minimize_free_energy(N, kappa, dx):
    res = minimize(edge_free_energy, initial_guess(N, dx), args=(N, kappa, dx))
    fields = res.x.reshape(2,2*N)
    return fields[0], fields[1]

if __name__ == '__main__':
    L = 5
    dx = 0.025
    N = int(L/dx)
    kappa = 1/math.sqrt(2)
    psi, ay = minimize_free_energy(N, kappa, dx)
    print edge_free_energy(np.append(psi, ay), N, kappa, dx)
    B = np.diff(ay)/dx
    f_dens = free_energy_density(psi, ay, N, kappa, dx)
    x = np.linspace(-N*dx, N*dx, 2*N, endpoint=False)
     
    plt.plot(x, psi)
#    plt.plot(x, ay)
    plt.plot((x + dx/2)[:-1], B*math.sqrt(2))
    plt.plot((x + dx/2)[:-1], f_dens)
    plt.xlim(min(x), max(x))
    plt.show()    
