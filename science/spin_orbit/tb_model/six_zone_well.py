from __future__ import division
import math
import matplotlib.pyplot as plt
import numpy as np
from tblib import *


def put_well_slice(ham, px, py, semiconductor_data, n_start, n_atoms):
    ts, tpar, tperp, P, Ec, Ev = semiconductor_data 
    for i in xrange(n_start, n_start + n_atoms):
        ham[0:2,i,0:2,i] = (Ec + 2*ts*(3 - math.cos(px) - math.cos(py)))*np.identity(2)

        Ax = np.array([[-1/math.sqrt(2), 0, 1/math.sqrt(6), 0],
                       [0, -1/math.sqrt(6), 0, 1/math.sqrt(2)]])
        Ay = np.array([[-1j/math.sqrt(2), 0, -1j/math.sqrt(6), 0],
                       [0, -1j/math.sqrt(6), 0, 1j/math.sqrt(2)]])
        Az = np.array([[0, math.sqrt(2/3), 0, 0],
                       [0, 0, math.sqrt(2/3), 0]])

        ham[0:2,i,2:,i] = (math.sqrt(2/3)*P*Ax*math.sin(px) 
                           + math.sqrt(2/3)*P*Ay*math.sin(py))
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

        if i > 0:
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
    put_well_slice(ham, px, py, hgte_data, n_cdte, n_hgte)
    put_well_slice(ham, px, py, cdte_data, n_cdte + n_hgte, n_cdte)

    return ham.reshape(6*N, 6*N)


if __name__ == '__main__':
    cdte_data  = 0.3, 1, 0.5, 0.9, 0.3, -0.5
    hgte_data  = 0.3, 1, 0.5, 0.9, -0.3, 0
    n_cdte = 30
    n_hgte = 5

    ham = lambda px: well_hamiltonian(px, 0, cdte_data, hgte_data, n_cdte, n_hgte)

    p_lim = (-math.pi/4, math.pi/4)
    pxrange, energies = stripe_energies(ham, (), NX = 81, plim = p_lim)

    for level in energies:
        plt.plot(pxrange, level)
    plt.xlim(*p_lim)
    plt.ylim(-0.7,0.4)

    plt.show()
