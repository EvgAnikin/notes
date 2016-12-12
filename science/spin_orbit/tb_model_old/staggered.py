from en_levels import Trans_matrix,energies
import numpy as np
import scipy.linalg as linalg
import math
import matplotlib.pyplot as plt

def staggered_hamiltonian(px,py,*args):
    ESO,xi,t1,t2,t3 = args

    H = np.zeros(2*3*2*3,dtype=complex).reshape((2,3,2,3))
    
    A,B = 0,1
    H[A,:,A,:] = np.diag([0,0,-ESO])
    H[B,:,B,:] = xi*np.diag(np.ones(3)) + np.diag([0,0,-ESO])
    H[A,:,B,:] = Trans_matrix(math.pi/4,t1,t2,t3)*(1+np.exp(-1j*(px+py))) +\
                   Trans_matrix(-math.pi/4,t1,t2,t3)*(np.exp(-1j*px) + np.exp(-1j*py))
    H[B,:,A,:] = Trans_matrix(math.pi/4,t1,t2,t3)*(1+np.exp(1j*(px+py))) +\
                   Trans_matrix(-math.pi/4,t1,t2,t3)*(np.exp(1j*px) + np.exp(1j*py))

    return H.reshape((2*3,2*3))

def slice_of_energies(ham, px, args, NX = 60):
    energies = []
    pyrange = np.linspace(0,2*math.pi, NX)
    
    for py in pyrange:
    	energies.append(linalg.eigvalsh(ham(px, py,*args)))
    						
    energies = np.array(energies).transpose()
    return pyrange, energies

def plot_slice(px, args):
    """ESO,xi,t1,t2,t3 = args"""
    PY, energies = slice_of_energies(staggered_hamiltonian,px,args)

    for level in energies:
        plt.plot(PY,level)
    plt.xlim(0,2*math.pi)
    plt.show()

def plot(args):
    """ESO,xi,t1,t2,t3 = args"""
    pxgrid, pygrid, es = energies(staggered_hamiltonian, args, NX = 40, NY = 40)
    	
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.set_xlim([-math.pi,math.pi])
    ax.set_ylim([-math.pi,math.pi])
        
    cmap = plt.get_cmap()
    for i in xrange(6):
        ax.plot_wireframe(pxgrid, pygrid, es[i], color = cmap(i/6.))
    
    plt.show()

if __name__ == '__main__':
    args = 1,0,0.4,0.2,0.15
    plot_slice(0.0,args)
