import math
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eigh
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


def hamiltonian(NX, NY, xi,m,t):
    H = np.zeros((NX*NY*2)**2, dtype=complex).reshape(NX,NY,2,NX,NY,2)  
    
    A,B = 0,1
    for i in xrange(NX):
        for j in xrange(NY):
            H[i,j,A,i,j,A] = xi
            H[i,j,B,i,j,B] = -xi

            if i < NX - 1:
                H[i,j,A,i+1,j,A] = 1./m/2
                H[i+1,j,A,i,j,A]   = 1./m/2
                H[i,j,B,i+1,j,B] = -1./m/2
                H[i+1,j,B,i,j,B] = -1./m/2

                H[i,j,A,i+1,j,B] = -1j*t
                H[i+1,j,A,i,j,B] = 1j*t
                H[i,j,B,i+1,j,A] = -1j*t
                H[i+1,j,B,i,j,A] = 1j*t
            
            if j < NY - 1:
                H[i,j,A,i,j+1,A] = 1./m/2
                H[i,j+1,A,i,j,A] = 1./m/2
                H[i,j,B,i,j+1,B] = -1./m/2
                H[i,j+1,B,i,j,B] = -1./m/2
    
                H[i,j,A,i,j+1,B] = 1*t
                H[i,j+1,A,i,j,B] = -1*t
                H[i,j,B,i,j+1,A] = -1*t
                H[i,j+1,B,i,j,A] = 1*t

    return H.reshape(NX*NY*2, NX*NY*2)


def closest_value_to_zero(energies):
    for i, e in enumerate(energies):
        if e >= 0:
           return i


def draw_state(vector, filename=None, show=False, magnitude_factor=5):
    rect_size = 10
    im_size = (rect_size*vector.shape[0], rect_size*vector.shape[1])
    im = Image.new('RGBA', im_size, (0,0,0,0))
    draw = ImageDraw.Draw(im)
    NX = vector.shape[0]
    max_val = np.amax(abs(vector[NX/2,0]))

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


if __name__ == '__main__':
    NX = 10
    NY = 10
    ham = hamiltonian(NX, NY, -0.4, 0.5, 1)
    energies,states = eigh(ham)
    
    states = np.transpose(states.reshape(NX,NY,2,NX*NY*2), (3,0,1,2))
    densities = np.sum(abs(states)**2, axis=3)
    
#    outfile = open('dens_small.npz', 'w')
#    np.savez(outfile, densities=densities)
    
    draw_state(densities[NX*NY], show=True)
    #cvz = closest_value_to_zero(energies)
    #print cvz
    #plot_state(0, vectors)
