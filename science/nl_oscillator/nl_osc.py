#encoding:utf-8
from __future__ import print_function
import math
import numpy as np
from numpy.linalg import eig, eigh, lstsq
import matplotlib.pyplot as plt


def nonlinear_oscillator(omega, beta, E, N):
    ham = np.zeros(N**2).reshape(N,N)

    for i in xrange(N):
        ham[i,i] = omega*i + beta/2*(i**2 - 1)
        if i < N - 1:
            ham[i,i+1] = math.sqrt(i+1)*E
            ham[i+1,i] = math.sqrt(i+1)*E

    return ham


def find_transition_matrix(omega_0, Omega, beta, E, N):
    ham = nonlinear_oscillator(omega_0 - Omega, beta, E, N)

    energies, states = eigh(ham)

    a_op = np.zeros(N**2).reshape(N,N)
    for i in xrange(N-1):
        a_op[i,i+1] = math.sqrt(i+1)

    a_matrix = np.zeros(N**2).reshape(N,N)
    for i in xrange(N):
        for j in xrange(N):
            a_matrix[i,j] = np.dot(states[:, i], np.dot(a_op, states[:,j]))

    P = np.zeros(N**2).reshape(N,N)
    for i in xrange(N):
        for j in xrange(N):
           gamma_ij = max(0, energies[j] - energies[i] + Omega)
#           gamma_ij = 0 if  energies[j] - energies[i] + Omega < 0 else 1
           P[i,j] = abs(gamma_ij*a_matrix[i,j])**2

    T = P - np.diag(np.sum(P, axis=0))
    return T, states


N = 40
T, states = find_transition_matrix(omega_0=1,
                           Omega=1.025,
                           beta=0.003,
                           E=0.0138,
                           N=N)

eigvals, vectors = eig(T)
numbered_eigvals = sorted(enumerate(eigvals), key = lambda tup: -tup[1])
equilibrium_vector = vectors[:, numbered_eigvals[0][0]] 

probabilities = [sum(np.abs(states[n,:])**2*np.abs(equilibrium_vector)) for n in xrange(N)]

#plt.plot(np.abs(vectors[:,numbered_eigvals[0][0]])**2)
#plt.plot(np.abs(vectors[:,numbered_eigvals[1][0]])**2)
ax = plt.subplot(111)
#ax.set_yscale('log')
plt.plot(probabilities)
plt.xlim(0,20)
plt.xlabel('Occupation number')
plt.ylabel('Probability')
#plt.show()
plt.savefig('bistability.png')
