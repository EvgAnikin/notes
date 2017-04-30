import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
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


def ham_with_spin(NX, NY, xi, m, t, imp_rate=0):
    H = np.zeros((NX*NY*2*2)**2, dtype=complex).reshape(NX,NY,2,2,NX,NY,2,2)  
    spin_up = 0
    spin_down = 1

    H[:, :, :, spin_up, :, :, :, spin_up] = hamiltonian(NX, NY, xi, m, t, 
                                                     spin=1, reshape=False)
    H[:, :, :, spin_down, :, :, :, spin_down] = hamiltonian(NX, NY, xi, m, t,
                                                         spin=-1, reshape=False)

    add_random_magnetic_impurities(H, imp_rate)
    return H.reshape(NX*NY*2*2, NX*NY*2*2)


def closest_value_to_zero(energies):
    for i, e in enumerate(energies):
        if e >= 0:
           return i


def draw_state(vector, filename=None, show=True, magnitude_factor=5):
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
            brightness = abs(vector[i,j]/max_val)
            color = (int(brightness*255), 0, 0)
            draw.rectangle([x, y, x + rect_size, y + rect_size], fill=color)
    
    im.save('new_fig.png' if not filename else filename)
    if show:
        im.show()


<<<<<<< HEAD
if __name__ == '__main__':
    NX = 10
    NY = 10
    ham = hamiltonian(NX, NY, -0.4, 0.5, 1)
=======
def compute_without_spin():
    NX = 10
    NY = 10
    ham = hamiltonian(NX, NY, -0.3, 0.4, 1)
>>>>>>> d07ae59250903a55249ee855c71de89d751ba3bc
    energies,states = eigh(ham)
    
    states = np.transpose(states.reshape(NX,NY,2,NX*NY*2), (3,0,1,2))
    densities = np.sum(abs(states)**2, axis=3)
    
    draw_state(densities[NX*NY], show=True)


#def compute_with_spin():
#    pass


if __name__ == '__main__':
    NX = 60
    NY = 8
    ham = ham_with_spin(NX, NY, -0.3, 1, 0.4, imp_rate=0.25)
    energies,states = eigh(ham)
    
    states = np.transpose(states.reshape(NX,NY,2,2,NX*NY*2*2), (4,0,1,2,3))
    densities = np.sum(abs(states)**2, axis=(3,4))
    
    print 'Number of states: {}'.format(NX*NY*2*2)
    draw_state(densities[2*NX*NY])
#    outfile = open('new_dens.npz', 'w')
#    np.savez(outfile, densities=densities)
