import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image, ImageDraw


def top_ins_twisted_bc(NX, NY, xi, m, t, p, spin=1, reshape=True):
    H = np.zeros((NX*NY*2)**2, dtype=complex).reshape(NX,NY,2,NX,NY,2)  
    
    A,B = 0,1
    px, py = p
    for i in xrange(NX):
        for j in xrange(NY):
            H[i,j,A,i,j,A] += xi + 2./m
            H[i,j,B,i,j,B] += -xi - 2./m

            x_twister = cmath.exp(1j*px) if i == NX - 1 else 1
              
            H[i,j,A,(i+1)%NX,j,A] += -1./m/2*x_twister
            H[(i+1)%NX,j,A,i,j,A] += -1./m/2/x_twister
            H[i,j,B,(i+1)%NX,j,B] += 1./m/2*x_twister
            H[(i+1)%NX,j,B,i,j,B] += 1./m/2/x_twister

            H[i,j,A,(i+1)%NX,j,B] += -1j*t*spin*x_twister
            H[(i+1)%NX,j,B,i,j,A] +=  1j*t*spin/x_twister
            H[i,j,B,(i+1)%NX,j,A] += -1j*t*spin*x_twister
            H[(i+1)%NX,j,A,i,j,B] +=  1j*t*spin/x_twister
            
            y_twister = cmath.exp(1j*py) if j == NY - 1 else 1

            H[i,j,A,i,(j+1)%NY,A] += -1./m/2*y_twister
            H[i,(j+1)%NY,A,i,j,A] += -1./m/2/y_twister
            H[i,j,B,i,(j+1)%NY,B] +=  1./m/2*y_twister
            H[i,(j+1)%NY,B,i,j,B] +=  1./m/2/y_twister
 
            H[i,j,A,i,(j+1)%NY,B] +=  1*t*y_twister
            H[i,(j+1)%NY,B,i,j,A] +=  1*t/y_twister
            H[i,j,B,i,(j+1)%NY,A] += -1*t*y_twister
            H[i,(j+1)%NY,A,i,j,B] += -1*t/y_twister
    
    return H.reshape(NX*NY*2, NX*NY*2) if reshape else H


def get_even_potential_disorder(NX, NY, magnitude):
    V = np.zeros((NX*NY*2)**2).reshape(NX,NY,2,NX,NY,2)

    for i in xrange(NX):
        for j in xrange(NY):
            potential = magnitude*(-1 + 2*np.random.rand())
            V[i, j, :, i, j, :] = potential*np.diag(np.ones(2))
            V[NX-i-1, NY-j-1, :, NX-i-1, NY-j-1, :] = potential*np.diag(np.ones(2))
    return V


def parity_top_ins(state):
    state_mod = np.copy(state)[::-1, ::-1, :]
    state_mod[:,:,1] *= -1

    epsilon = 1e-6
    if np.sum(abs(state - state_mod)**2) < epsilon:
        return 1
    elif np.sum(abs(state + state_mod)**2) < epsilon:
        return -1
    else:
        return None


def get_topological_invariant():
    NX, NY, xi, m, t = 15, 15, 0.1, 1, 0.4

    time_reversal_invariant_momenta = [(0,0), (0,math.pi), (math.pi,0), (math.pi, math.pi)]
    invariant = 1
    V = get_even_potential_disorder(NX,NY,1)
    energies_for_histogram = []
    for px, py in time_reversal_invariant_momenta:
        ham = top_ins_twisted_bc(NX, NY, xi, m, t, (px, py), reshape=False) + V
        ham = ham.reshape(NX*NY*2, NX*NY*2)
        #print ham
        energies, states = eigh(ham)
        energies_for_histogram.append(energies)
        states = states.T.reshape(NX*NY*2, NX, NY, 2)
        for e, s in zip(energies[:NX*NY], states[:NX*NY]):
            p = parity_top_ins(s)
            if not p:
                return None
            else:
                invariant *= p
    return invariant


if __name__ == '__main__':
    print get_topological_invariant()
