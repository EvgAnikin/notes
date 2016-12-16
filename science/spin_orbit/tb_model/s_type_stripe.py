import math
import numpy as np
import matplotlib.pyplot as plt
from tblib import *


def staggered_stripe_hamiltonian(py,N,ESO,xi,t1,t2,t3):

    H = np.zeros((2*3*N)**2, dtype=complex).reshape(N,2,3,N,2,3)
    A,B = 0,1
    for i in xrange(N):
        H[i,A,:,i,A,:] = np.diag([0,0,-ESO])
        H[i,B,:,i,B,:] = np.diag([0,0,-ESO]) + xi*np.diag(np.ones(3))
        H[i,A,:,i,B,:] = Trans_matrix(math.pi/4,t1,t2,t3) + \
                            np.exp(-1j*py)*Trans_matrix(-math.pi/4,t1,t2,t3)
        H[i,B,:,i,A,:] = Trans_matrix(math.pi/4,t1,t2,t3) + \
                            np.exp(1j*py)*Trans_matrix(-math.pi/4,t1,t2,t3)
        if i > 0:
            H[i,A,:,i-1,B,:] = np.exp(-1j*py)*Trans_matrix(math.pi/4,t1,t2,t3) + \
                            Trans_matrix(-math.pi/4,t1,t2,t3)
            H[i-1,B,:,i,A,:] = np.exp(1j*py)*Trans_matrix(math.pi/4,t1,t2,t3) + \
                            Trans_matrix(-math.pi/4,t1,t2,t3)

    return H.reshape(2*3*N,2*3*N)

def s_type_stripe_hamiltonian(py,N,ts,tsp,Es,tlong, ttrans,t3,ESO):
    H = np.zeros((4*N)**2, dtype = complex).reshape(N,4,N,4)

    for i in xrange(N):
        H[i,:,i,:] += np.diag([Es + 2*ts*(2 - np.cos(py)),\
                    2*ttrans*np.cos(py),\
                    2*tlong*np.cos(py),\
                    2*t3*np.cos(py)])
        H[i,0,i,2] += 2*tsp*(-1j)*np.sin(py)
        H[i,2,i,0] += 2*tsp*1j*np.sin(py)

        H[i,1:,i,1:] += -ESO/3. * np.array([(1,1j,-1), \
						(-1j,1,1j), \
						(-1, -1j, 1.)])

        if i > 0:
            H[i,:,i-1,:] += np.array( [ (-ts,tsp,0,0),\
                                       (-tsp,tlong,0,0),\
                                       (0, 0,ttrans,0),\
                                       (0, 0,0, t3) ])
            H[i-1,:,i,:] += np.array( [ (-ts,-tsp,0,0),\
                                       (tsp,tlong,0,0),\
                                       (0, 0,ttrans,0),
                                       (0, 0,0, t3) ])
    
    H = H.reshape(4*N,4*N)
    return H


def plot(args):
    """ args = N,ts,tsp,Es,tlong, ttrans,t3,ESO,td"""
    PY, energies = stripe_energies(s_type_stripe_hamiltonian,args,\
                                    plim = (-math.pi,math.pi),NX = 80)

    for level in energies:
        plt.plot(PY,level)
    plt.xlim(-math.pi,math.pi)
    plt.show()

if __name__ == '__main__':
#    args = 30,0.2,0.15,1.6,0.3,0.15,0.10,4
#    args = 30,0.2,0.2,1.,0.5,0.3,0.15,1
    args = 30,0.2,0.2,1,0.4,0.2,0.15,1
    
#    args = 30,0.2,0,1,0,0,0,0
    plot(args)
