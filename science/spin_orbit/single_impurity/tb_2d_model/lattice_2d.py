import math
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eigh
from mpl_toolkits.mplot3d import Axes3D


def one_band_ham(NX, NY, t):
    H = np.zeros((NX*NY)**2, dtype=complex).reshape(NX,NY,NX,NY)  
    
    for i in xrange(NX):
        for j in xrange(NY):
            H[i,j,(i+1)%NX,j] = -t
            H[i,j,(i-1)%NX,j] = -t
            H[i,j,i,(j+1)%NY] = -t
            H[i,j,i,(j-1)%NY] = -t

    H[0,0,0,0] = -5.

    return H.reshape(NX*NY, NX*NY)


def hamiltonian(NX, NY, xi,m,t):
    H = np.zeros((NX*NY*2)**2, dtype=complex).reshape(NX,NY,2,NX,NY,2)  
    
    A,B = 0,1
    for i in xrange(NX):
        for j in xrange(NY):
            H[i,j,A,i,j,A] = xi
            H[i,j,B,i,j,B] = -xi

            H[i,j,A,(i+1)%NX,j,A] = 1./m/2
            H[i,j,A,(i-1)%NX,j,A] = 1./m/2
            H[i,j,A,i,(j+1)%NY,A] = 1./m/2
            H[i,j,A,i,(j-1)%NY,A] = 1./m/2

            H[i,j,B,(i+1)%NX,j,B] = -1./m/2
            H[i,j,B,(i-1)%NX,j,B] = -1./m/2
            H[i,j,B,i,(j+1)%NY,B] = -1./m/2
            H[i,j,B,i,(j-1)%NY,B] = -1./m/2

            H[i,j,A,(i+1)%NX,j,B] = -1j
            H[i,j,A,(i-1)%NX,j,B] = 1j
            H[i,j,A,i,(j+1)%NY,B] = 1
            H[i,j,A,i,(j-1)%NY,B] = -1

            H[i,j,B,(i+1)%NX,j,A] = -1j
            H[i,j,B,(i-1)%NX,j,A] = 1j
            H[i,j,B,i,(j+1)%NY,A] = -1
            H[i,j,B,i,(j-1)%NY,A] = 1

    return H.reshape(NX*NY*2, NX*NY*2)


def closest_value_to_zero(energies):
    for i, e in enumerate(energies):
        if e >= 0:
           return i


def plot_state(number):
    state = vectors[:,number]
    state_mod = state.reshape(NX,NY)
    
    plt.plot(range(0,NX), state_mod[:,0])#map(abs, state_mod[:,0,0]))
    plt.show()


NX = 10
NY = 15
ham = one_band_ham(NX,NY,1)
energies,vectors = eigh(ham)

#cvz = closest_value_to_zero(energies)
#print cvz
plot_state(0)
