from __future__ import division
import math
import matplotlib.pyplot as plt
import numpy as np
from tblib import *


def well_hamiltonian(px, py, ts, tpar, tperp, P, Ec, Ev, N):
    
    ham = np.zeros((6*N)**2, dtype=complex).reshape(6,N,6,N)

    for i in xrange(N):
        ham[0:2,i,0:2,i] = (Ec + 2*ts*(3 - math.cos(px) - math.cos(py)))*np.identity(2)

        Ax = 1./2*np.array([[-math.sqrt(3),0,1,0],
                             [0,-1,0,math.sqrt(3)]],dtype=complex)
        Ay = -1.j/2*np.array([[math.sqrt(3),0,1,0],
                             [0,1,0,math.sqrt(3)]],dtype=complex)
        Az = np.array([[0,1,0,0],
                       [0,0,1,0]], dtype=complex)
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
    
            ham[0,i-1,3,i] = 1j/2*math.sqrt(2/3)*P
            ham[3,i,0,i-1] = -1j/2*math.sqrt(2/3)*P
    
            ham[1,i-1,4,i] = 1j/2*math.sqrt(2/3)*P
            ham[4,i,1,i-1] = -1j/2*math.sqrt(2/3)*P
    
            ham[2,i-1,2,i] = ham[2,i,2,i-1] = ham[5,i-1,5,i] = ham[5,i,5,i-1] = tperp
            ham[3,i-1,3,i] = ham[3,i,3,i-1] = ham[4,i-1,4,i] = ham[4,i,4,i-1] =\
                                            (tperp + 2*tpar)/3

    return ham.reshape(6*N, 6*N)


if __name__ == '__main__':
    ts, tpar, tperp, P, Ec, Ev = 0.3, 1, 0.5, 0.9, -0.5, 0
    N = 10
    ham = lambda px: well_hamiltonian(px, 0, ts, tpar, tperp, P, Ec, Ev, N)

    pxrange, energies = stripe_energies(ham, (), NX = 61, plim = (-math.pi, math.pi))

    for level in energies:
        plt.plot(pxrange, level)
    plt.xlim(-math.pi, math.pi)

    plt.show()