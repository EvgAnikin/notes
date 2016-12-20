import math
import matplotlib.pyplot as plt
import numpy as np
from tblib import stripe_energies

h2m = 0.038

def p_band_ham(Ev, gamma1, gamma2, kx, ky, kz):
    Jx = np.array([[0,math.sqrt(3)/2,0,0],
                  [math.sqrt(3)/2,0,1,0],
                  [0,1,0,math.sqrt(3)/2],
                  [0,0,math.sqrt(3)/2,0]], dtype=complex)
    Jy = np.array([[0,-1j*math.sqrt(3)/2,0,0],
                  [1j*math.sqrt(3)/2,0,-1j,0],
                  [0,1j,0,-1j*math.sqrt(3)/2],
                  [0,0,1j*math.sqrt(3)/2,0]], dtype=complex)
    Jz = np.diag([3./2 + 0j,1./2,-1./2,-3./2])

#    Vx = np.dot(Jx, np.dot(Jy,Jy) - np.dot(Jz,Jz))
#    Vy = np.dot(Jy, np.dot(Jz,Jz) - np.dot(Jx,Jx))
#    Vz = np.dot(Jz, np.dot(Jx,Jx) - np.dot(Jy,Jy))

    return (Ev*np.identity(4) + 
            h2m*(-(gamma1 + 5./2*gamma2)*(kx**2 + ky**2 + kz**2)*np.identity(4) +
            2*gamma2*np.dot((kx*Jx + ky*Jy + kz*Jz), (kx*Jx + ky*Jy + kz*Jz))))
#            4*(gamma3-gamma2)*
#                (kx**2*np.dot(Jx, Jx) + ky**2*np.dot(Jy, Jy) + kz**2*np.dot(Jz, Jz)))
#            t*(kx*Vx + ky*Vy + kz*Vz))

def kane_ham(kx, ky, kz, args):
    Es, Ev, F, P, gamma1, gamma2 = args
    Ax = 1./2*np.array([[-math.sqrt(3),0,1,0],
                         [0,-1,0,math.sqrt(3)]],dtype=complex)
    Ay = -1.j/2*np.array([[math.sqrt(3),0,1,0],
                         [0,1,0,math.sqrt(3)]],dtype=complex)
    Az = np.array([[0,1,0,0],
                   [0,0,1,0]], dtype=complex)

    H = np.zeros(6*6, dtype=complex).reshape(6,6)
    H[0:2,0:2] = np.identity(2)*(Es + h2m*(2*F + 1)*(kx**2 + ky**2 + kz**2))
    H[2:,2:] = p_band_ham(Ev, gamma1, gamma2, kx, ky, kz)
    H[0:2,2:] = math.sqrt(2./3)*P*(kx*Ax + ky*Ay + kz*Az)
    H[2:,0:2] = H[0:2,2:].conj().T

    return H


ham_x = lambda pz: kane_ham(0, 0, pz, (1.606-0.57, -0.57, -0.09, 0.845, 1.47, -0.28))

pyrange, energies = stripe_energies(ham_x, (), NX = 101, plim = (-1,1))

for band in energies:
    plt.plot(pyrange, band)

plt.ylim(-1.5,1.5)
plt.show()
