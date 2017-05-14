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


def one_band_ham(NX, NY, t):
    H = np.zeros((NX*NY)**2, dtype=complex).reshape(NX,NY,NX,NY)  
    
    for i in xrange(NX):
        for j in xrange(NY):
            if i < NX - 1:
                H[i,j,i+1,j] = -t
                H[i+1,j,i,j] = -t
            if j < NY - 1:
                H[i,j,i,j+1] = -t
                H[i,j+1,i,j] = -t

    return H.reshape(NX*NY, NX*NY)


def hamiltonian(NX, NY, xi, m, t, spin=1, reshape=True):
    H = np.zeros((NX*NY*2)**2, dtype=complex).reshape(NX,NY,2,NX,NY,2)  
    
    A,B = 0,1
    for i in xrange(NX):
        for j in xrange(NY):
            H[i,j,A,i,j,A] = xi
            H[i,j,B,i,j,B] = -xi

#            if i < NX - 1:
            H[i,j,A,(i+1)%NX,j,A] = 1./m/2
            H[(i+1)%NX,j,A,i,j,A]   = 1./m/2
            H[i,j,B,(i+1)%NX,j,B] = -1./m/2
            H[(i+1)%NX,j,B,i,j,B] = -1./m/2

            H[i,j,A,(i+1)%NX,j,B] = -1j*t*spin
            H[(i+1)%NX,j,A,i,j,B] = 1j*t*spin
            H[i,j,B,(i+1)%NX,j,A] = -1j*t*spin
            H[(i+1)%NX,j,B,i,j,A] = 1j*t*spin
            
            if j < NY - 1:
                H[i,j,A,i,j+1,A] = 1./m/2
                H[i,j+1,A,i,j,A] = 1./m/2
                H[i,j,B,i,j+1,B] = -1./m/2
                H[i,j+1,B,i,j,B] = -1./m/2
    
                H[i,j,A,i,j+1,B] = 1*t
                H[i,j+1,A,i,j,B] = -1*t
                H[i,j,B,i,j+1,A] = -1*t
                H[i,j+1,B,i,j,A] = 1*t
    
    return H.reshape(NX*NY*2, NX*NY*2) if reshape else H


def ham_with_obstacle_1(NX, NY, xi, m, t, spin=1, reshape=True):
    H = hamiltonian(NX, NY, xi, m, t, reshape=False)
    
    length = 1
    A, B = 0, 1
    for i in xrange(1, length+1):
        H[NX/2, i, A, NX/2, i, A] = H[NX/2, i, B, NX/2, i, B] = 100/m
    return H.reshape(NX*NY*2, NX*NY*2)


def ham_with_obstacle_2(NX, NY, xi, m, t, spin=1, reshape=True):
    H = hamiltonian(NX, NY, xi, m, t, reshape=False)
    
    length = 2
    A, B = 0, 1
    for i in xrange(0, length):
        x = NX/2 - length/2 + i
        y = 0
        H[x, y, A, x, y, A] = H[x, y, B, x, y, B] = 100/m
    return H.reshape(NX*NY*2, NX*NY*2)
    


def get_random_magnetic_impurity():
    phi = 2*math.pi*np.random.rand()
    theta = math.acos(-1 + 2*np.random.rand())

    imp = np.zeros(2**4, dtype=complex).reshape(2,2,2,2)
    imp[0,:,0,:] = imp[1,:,1,:] = np.array([[math.cos(theta), 
                                             math.sin(theta)*cmath.exp(-1j*phi)],
                                            [math.sin(theta)*cmath.exp(1j*phi),
                                            -math.cos(theta)]])
    return imp
    

def add_random_magnetic_impurities(H, magnitude):
    NX = H.shape[0]
    NY = H.shape[1]

    for i in xrange(NX):
        for j in xrange(NY):
            H[i, j, :, :, i, j, :, :] += magnitude*get_random_magnetic_impurity()


def add_potential_disorder(H, magnitude):
    NX = H.shape[0]
    NY = H.shape[1]

    for i in xrange(NX):
        for j in xrange(NY):
            potential = magnitude*(-1 + 2*np.random.rand())
            H[i, j, :, :, i, j, :, :] += potential*np.diag(np.ones(4)).reshape(2,2,2,2)


def ham_with_spin(NX, NY, xi, m, t, pot_imp_rate=0, mag_imp_rate=0):
    H = np.zeros((NX*NY*2*2)**2, dtype=complex).reshape(NX,NY,2,2,NX,NY,2,2)  
    spin_up = 0
    spin_down = 1

    H[:, :, :, spin_up, :, :, :, spin_up] = hamiltonian(NX, NY, xi, m, t, 
                                                     spin=1, reshape=False)
    H[:, :, :, spin_down, :, :, :, spin_down] = hamiltonian(NX, NY, xi, m, t,
                                                         spin=-1, reshape=False)

    add_potential_disorder(H, pot_imp_rate)
    add_random_magnetic_impurities(H, mag_imp_rate)
    return H.reshape(NX*NY*2*2, NX*NY*2*2)


def closest_value_to_zero(energies):
    for i, e in enumerate(energies):
        if e >= 0:
           return i


def draw_state(vector, filename=None, show=True, magnitude_factor=3):
    rect_size = 10
    im_size = (rect_size*vector.shape[0], rect_size*vector.shape[1])
    im = Image.new('RGBA', im_size, (0,0,0,0))
    draw = ImageDraw.Draw(im)
    NX = vector.shape[0]
    max_val = np.amax(abs(vector))

    nx, ny = vector.shape
    for i in xrange(nx):
        for j in xrange(ny):
            x = i*rect_size
            y = j*rect_size
        
#            brightness = max(0, 1 + math.log(abs(vector[i,j]/max_val), 10)/magnitude_factor)
            brightness = min(1, abs(vector[i,j]/max_val))
            color = (int(brightness*255), 0, 0)
            draw.rectangle([x, y, x + rect_size, y + rect_size], fill=color)
    
    im.save('new_fig.png' if not filename else filename)
    if show:
        im.show()


if __name__ == '__main__':
    NX = 10
    NY = 10
    ham = hamiltonian(NX, NY, -0.4, 0.5, 1)


def compute_without_spin():
    NX = 10
    NY = 10
    ham = hamiltonian(NX, NY, -0.3, 0.4, 1)
    energies,states = eigh(ham)
    
    states = np.transpose(states.reshape(NX,NY,2,NX*NY*2), (3,0,1,2))
    densities = np.sum(abs(states)**2, axis=3)
    
    draw_state(densities[NX*NY], show=True)


def plot_histogram(energies, bins):
    values, base = np.histogram(energies, bins=bins)
    plt.plot((base[:-1] + base[1:])/2, values)
    plt.ylim = (-2,2)
    plt.show()


def diag_with_obstacle():
    NX = 16
    NY = 16
    ham = ham_with_obstacle_1(NX, NY, -0.2, 1, 0.4)
    
    n_of_states = NX*NY*2
    lo, hi = NX*NY - n_of_states/2, NX*NY + n_of_states/2 - 1
    t0 = time.clock()
    energies,states = eigh(ham, eigvals=(lo, hi))
    print 'time: {}'.format(time.clock() - t0)
    print 'NX: {} NY: {}'.format(NX, NY)
#    plot_histogram(energies, bins=20)
    
    states = np.transpose(states.reshape(NX,NY,2,n_of_states), (3,0,1,2))
    densities = np.sum(abs(states)**2, axis=3)
    return densities


def diagonalize_w_spin_and_imp():
    NX = 10
    NY = 20
    ham = ham_with_spin(NX, NY, -0.2, 1, 0.4, pot_imp_rate=0.0, mag_imp_rate=0.25)
    
    n_of_states = NX*NY*2*2
    lo, hi = NX*NY*2 - n_of_states/2, NX*NY*2 + n_of_states/2 - 1
    t0 = time.clock()
    energies,states = eigh(ham, eigvals=(lo, hi))
    print 'time: {}'.format(time.clock() - t0)
    plot_histogram(energies, bins=20)
    
    states = np.transpose(states.reshape(NX,NY,2,2,hi-lo), (4,0,1,2,3))
    densities = np.sum(abs(states)**2, axis=(3,4))
    return densities


if __name__ == '__main__':
    densities = diag_with_obstacle()

#def diag_with_obstacle():
    
    
#    print 'Number of states: {}'.format(NX*NY*2*2)
#    draw_state(densities[2*NX*NY])
#    outfile = open('new_dens.npz', 'w')
#    np.savez(outfile, densities=densities)
