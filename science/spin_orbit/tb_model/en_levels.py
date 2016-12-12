#encoding:utf-8
import math
import numpy as np
import numpy.linalg as linalg
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from tblib import *

def hamiltonian(px,py, tlong, ttrans,t3,ESO,td):
    "hamiltonian(px,py, tlong, ttrans,t3,ESO,td)"
    
    H0 = 2*np.array(\
    			[ (tlong*np.cos(px) + ttrans*np.cos(py),\
    				4*td*np.sin(px)*np.sin(py), 0),\
    			  (4*td*np.sin(px)*np.sin(py),\
    				 tlong*np.cos(py) + ttrans*np.cos(px), 0),\
    			  (0,0,t3*(np.cos(px) + np.cos(py))) ]\
    			)
    
    V =	-ESO/3. * np.array([(1,1j,-1), \
						(-1j,1,1j), \
						(-1, -1j, 1.)])
	
    return H0 + V

def s_type_hamiltonian(px,py,ts,tsp,Es,tlong,ttrans,t3,ESO):
    "hamiltonian(px,py,ts,tlong,ttrans,t3,ESO)"
    
    H0 = np.zeros(16,dtype = complex).reshape((4,4))
    H0 += np.diag([Es + 2*ts*(2 - np.cos(px) - np.cos(py)),\
                    2*(tlong*np.cos(px) + ttrans*np.cos(py)),\
                    2*(tlong*np.cos(py) + ttrans*np.cos(px)),\
                    2*t3*(np.cos(px) + np.cos(py))])

    H0[0,1:3] += 2*tsp*(-1j)*np.array([np.sin(px),np.sin(py)])
    H0[1:3,0] += 2*tsp*1j*np.array([np.sin(px),np.sin(py)])
    
    V =	-ESO/3. * np.array([(0,0,0,0),\
                        (0,1,1j,-1), \
						(0,-1j,1,1j), \
						(0,-1, -1j, 1.)])
	
    return H0 + V

def hamiltonian_2(px,py, tlong, ttrans,t3,ESO,td):
	"hamiltonian_2(px,py, tlong, ttrans,t3,ESO,td)"
	t1 = 1./2*(tlong + ttrans)
	t2 = 1./2*(tlong - ttrans)

	atomic = np.diag([0,0,-ESO])	

	return atomic + 2*np.cos(px)*Trans_matrix(0,t1,t2,t3) + 2*np.cos(py)*Trans_matrix(np.pi/2,t1,t2,t3)

def hamiltonian_3(px,py, ESO,t1,t2,t3,td1,td2,td3):
	"hamiltonian_3(px,py,ESO,t1,t2,t3,td1,td2,td3)"
	atomic = np.diag([0,0,-ESO])	

	return atomic + 2*np.cos(px)*Trans_matrix(0,t1,t2,t3) + 2*np.cos(py)*Trans_matrix(np.pi/2,t1,t2,t3) +\
					2*np.cos(px+py)*Trans_matrix(np.pi/4,td1,td2,td3) + 2*np.cos(px-py)*Trans_matrix(np.pi/4,td1,td2,td3)


if __name__ == '__main__':
    args = 0.2,0.2,1.,0.5,0.3,0.15,1
    pxgrid, pygrid, es = energies(s_type_hamiltonian, args, NX = 40, NY = 40)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.set_xlim([-math.pi,math.pi])
    ax.set_ylim([-math.pi,math.pi])
#    ax.plot_wireframe(pxgrid, pygrid, es[1]-es[0], color = 'red')
#    ax.contour(pxgrid, pygrid, es[2]-es[1], 20)
    
    ax.plot_wireframe(pxgrid, pygrid, es[0], color = 'blue')
    ax.plot_wireframe(pxgrid, pygrid, es[1], color = 'red')
    ax.plot_wireframe(pxgrid, pygrid, es[2], color = 'green')
    ax.plot_wireframe(pxgrid, pygrid, es[3], color = 'orange')
    
#    ax.contour(pxgrid, es[0], pygrid, 30)
#    ax.contour(pxgrid, es[1], pygrid, 30)
#    ax.contour(pxgrid, es[2], pygrid, 30)
    
    plt.show()
