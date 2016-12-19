# coding: utf-8
import math
import matplotlib.pyplot as plt
import numpy as np
from tblib import stripe_energies


def p_band_ham(a,b,c,t,kx,ky,kz):
    Jx = np.array([[0,math.sqrt(3)/2,0,0],
                  [math.sqrt(3)/2,0,1,0],
                  [0,1,0,math.sqrt(3)/2],
                  [0,0,math.sqrt(3)/2,0]], dtype=complex)
    Jy = np.array([[0,-1j*math.sqrt(3)/2,0,0],
                  [1j*math.sqrt(3)/2,0,-1j,0],
                  [0,1j,0,-1j*math.sqrt(3)/2],
                  [0,0,1j*math.sqrt(3)/2,0]], dtype=complex)
    Jz = np.diag([3./2 + 0j,1./2,-1./2,-3./2])

    Vx = np.dot(Jx, np.dot(Jy,Jy) - np.dot(Jz,Jz))
    Vy = np.dot(Jy, np.dot(Jz,Jz) - np.dot(Jx,Jx))
    Vz = np.dot(Jz, np.dot(Jx,Jx) - np.dot(Jy,Jy))

    return (a*(kx**2 + ky**2 + kz**2)*np.identity(4) +
            b*np.dot((kx*Jx + ky*Jy + kz*Jz),(kx*Jx + ky*Jy + kz*Jz)) +
            c*(kx**2*np.dot(Jx,Jx) + ky**2*np.dot(Jy,Jy) + kz**2*np.dot(Jz,Jz)) +
            t*(kx*Vx + ky*Vy + kz*Vz))

a,b,c,t = -1,1,0.5,0.0
ham_x = lambda px: ham(a,b,c,t,px,0.4,0)

pyrange, energies = stripe_energies(ham_x, (), NX = 101, plim = (-1,1))

for band in energies:
    plt.plot(pyrange, band)

plt.show()
