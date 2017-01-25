from __future__ import division
import math
import matplotlib.pyplot as plt
import numpy as np
from tblib import *


def put_well_slice(ham, px, py, semiconductor_data, n_start, n_atoms, is_middle=False):
    ts, tpar, tperp, P, Ec, Ev = semiconductor_data 

    Ax = np.array([[-1/math.sqrt(2), 0, 1/math.sqrt(6), 0],
                   [0, -1/math.sqrt(6), 0, 1/math.sqrt(2)]])
    Ay = np.array([[-1j/math.sqrt(2), 0, -1j/math.sqrt(6), 0],
                   [0, -1j/math.sqrt(6), 0, -1j/math.sqrt(2)]])
    Az = np.array([[0, math.sqrt(2/3), 0, 0],
                   [0, 0, math.sqrt(2/3), 0]])

    for i in xrange(n_start, n_start + n_atoms):
        ham[0:2,i,0:2,i] = (Ec + 2*ts*(3 - math.cos(px) - math.cos(py)))*np.identity(2)

        ham[0:2,i,2:,i] = P*(Ax*math.sin(px) + Ay*math.sin(py))
        ham[2:,i,0:2,i] = ham[0:2,i,2:,i].T.conj()
        
        ham[2:,i,2:,i] =  ((Ev - 2*tpar - 4*tperp)*np.identity(4) + 
                          np.array([[(tpar + tperp)*(math.cos(px) + math.cos(py)), 0,
                            -1/math.sqrt(3)*(tpar - tperp)*(math.cos(px) - math.cos(py)), 0],
                            [0, (tpar/3 + 5*tperp/3)*(math.cos(px) + math.cos(py)), 
                             0, -1/math.sqrt(3)*(tpar - tperp)*(math.cos(px) - math.cos(py))],
                            [-1/math.sqrt(3)*(tpar - tperp)*(math.cos(px) - math.cos(py)), 0,
                             (tpar/3 + 5*tperp/3)*(math.cos(px) + math.cos(py)), 0],
                            [0, -1/math.sqrt(3)*(tpar - tperp)*(math.cos(px) - math.cos(py)),
                             0, (tpar + tperp)*(math.cos(px) + math.cos(py))]]))

    if is_middle:
        range_of_hopping = (n_start, n_start + n_atoms + 1) 
    else:
        range_of_hopping = (n_start + 1, n_start + n_atoms) 

    for i in xrange(*range_of_hopping):
        ham[0,i-1,0,i] = ham[0,i,0,i-1] = -ts
        ham[1,i-1,1,i] = ham[1,i,1,i-1] = -ts
  
        ham[0:2,i-1,2:,i] = 1j*P/2*Az
        ham[0:2,i,2:,i-1] = -1j*P/2*Az
  
        ham[2:,i,0:2,i-1] = ham[0:2,i-1,2:,i].T.conj()
        ham[2:,i-1,0:2,i] = ham[0:2,i,2:,i-1].T.conj()
  
        ham[2,i-1,2,i] = ham[2,i,2,i-1] = ham[5,i-1,5,i] = ham[5,i,5,i-1] = tperp
        ham[3,i-1,3,i] = ham[3,i,3,i-1] = ham[4,i-1,4,i] = ham[4,i,4,i-1] =\
                                      (tperp + 2*tpar)/3


def well_hamiltonian(px, py, cdte_data, hgte_data, n_cdte, n_hgte):
    N = 2*n_cdte + n_hgte
    ham = np.zeros((6*N)**2, dtype=complex).reshape(6,N,6,N)

    put_well_slice(ham, px, py, cdte_data, 0, n_cdte)
    put_well_slice(ham, px, py, hgte_data, n_cdte, n_hgte, is_middle=True)
    put_well_slice(ham, px, py, cdte_data, n_cdte + n_hgte, n_cdte)

    return ham.reshape(6*N, 6*N)


def get_tb_params(gamma1, gamma2, Ep, Ec, Ev):
    hbar = 1.0545718e-27
    me = 9.10938356e-28
    eV = 1.6021766209e-12
    a = 6.46e-8
    E0 = hbar**2/2./me/a**2/eV

    tpar = E0*(gamma1 + 4*gamma2)
    tperp = E0*(gamma1 - 2*gamma2)
    ts = E0
    P = math.sqrt(Ep*E0)

    return ts, tpar, tperp, P, Ec, Ev


if __name__ == '__main__':
    cdte_data = get_tb_params(
                        gamma1=1.47,
                        gamma2=-0.28,
                        Ep=18.8,
                        Ec=1.606-0.57,
                        Ev=-0.57)
    hgte_data = get_tb_params(
                        gamma1=4.1,
                        gamma2=0.5,
                        Ep=18.8,
                        Ec=-0.303,
                        Ev=0)
    n_cdte = 40

    for width in xrange(1,15):
        ham = lambda px: well_hamiltonian(px, 0, cdte_data, hgte_data, 
                                          n_cdte=n_cdte, n_hgte=width)
        energies = linalg.eigvalsh(ham(0))
        plt.scatter(width*np.ones(len(energies)), energies, marker = '.')

    plt.ylim(-1,1)
    plt.show()
