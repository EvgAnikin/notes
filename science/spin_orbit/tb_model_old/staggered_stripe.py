from stripe import stripe_energies
from en_levels import Trans_matrix
import math
import numpy as np
import matplotlib.pyplot as plt

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

def plot(args):
    """ args = N,ESO,xi,t1,t2,t3 """
    PY, energies = stripe_energies(staggered_stripe_hamiltonian,args)

    for level in energies:
        plt.plot(PY,level)
    plt.xlim(0,2*math.pi)
    plt.show()

if __name__ == '__main__':
    args = 30,6,4,3,1,1
    plot(args)
