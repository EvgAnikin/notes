import math
import numpy as np
import numpy.linalg as linalg
import math

def Trans_matrix(phi, *args): 
	t1,t2,t3 = args
	return  np.array( [(t1, 1/np.sqrt(3)*t2*np.exp(-2j*phi),\
									 np.sqrt(2./3)*t2*np.exp(-2j*phi)),\
						(1/np.sqrt(3)*t2*np.exp(2j*phi), 1./3*(t1+2*t3), np.sqrt(2)/3*(t1-t3)),\
						(np.sqrt(2./3)*t2*np.exp(2j*phi),np.sqrt(2)/3*(t1-t3), 1./3*(2*t1+t3))] )

def energies(ham_2d, args, NX = 40, NY = 40):
    """Returns the array of energies associated with each of the wave vector from the
    2D Brillouin zone. The hamiltonian must be of form ham(px,py,*args)""" 
    energies = []
    pxrange = np.linspace(-math.pi,math.pi, NX)
    pyrange = np.linspace(-math.pi,math.pi, NY)
    
    for px in pxrange:
    	energies.append([])
    	for py in pyrange:
    		energies[-1].append(linalg.eigvalsh(ham_2d(px,py,*args)))
    
    pxgrid, pygrid = np.meshgrid(pxrange,pyrange,indexing = 'ij')
    						
    energies = np.array(energies).transpose(2,0,1)
    return pxgrid, pygrid, energies

def stripe_energies(stripe_ham, args, NX = 60, plim = (0,2*math.pi)):
    energies = []
    pyrange = np.linspace(plim[0],plim[1], NX)
    
    for py in pyrange:
    	energies.append(linalg.eigvalsh(stripe_ham(py,*args)))
    						
    energies = np.array(energies).transpose()
    return pyrange, energies

def slice_of_energies(ham, px, args, NX = 60):
    energies = []
    pyrange = np.linspace(0,2*math.pi, NX)
    
    for py in pyrange:
    	energies.append(linalg.eigvalsh(ham(px, py,*args)))
    						
    energies = np.array(energies).transpose()
    return pyrange, energies
