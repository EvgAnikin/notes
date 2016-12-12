#encoding:utf-8
import math
import numpy as np
import numpy.linalg as linalg
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

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


#тот же гамильтониан в другом базисе
#def hamiltonian_1(px,py, tlong, ttrans,t3,ESO,td):
#	t1 = 1./2*(tlong + ttrans)
#	t2 = 1./2*(tlong - ttrans)
#
#	atomic = np.diag([0,0,-ESO])	
#	tx_matrix = np.array([ (t1, t2/math.sqrt(3), math.sqrt(2./3)*t2),\
#						(t2/math.sqrt(3), (t1+2*t3)/3., math.sqrt(2)/3*(t1-t3)),\
#						(math.sqrt(2./3)*t2, math.sqrt(2)/3*(t1-t3), t1) ])
#	ty_matrix = np.array([ (t1, -t2/math.sqrt(3), -math.sqrt(2./3)*t2),\
#						(-t2/math.sqrt(3), (t1+2*t3)/3., math.sqrt(2)/3*(t1-t3)),\
#						(-math.sqrt(2./3)*t2, math.sqrt(2)/3*(t1-t3), t1) ])
#	return atomic + 2*math.cos(px)*tx_matrix + 2*math.cos(py)*ty_matrix

def Trans_matrix(phi, *args): 
	t1,t2,t3 = args
	return  np.array( [(t1, 1/np.sqrt(3)*t2*np.exp(-2j*phi),\
									 np.sqrt(2./3)*t2*np.exp(-2j*phi)),\
						(1/np.sqrt(3)*t2*np.exp(2j*phi), 1./3*(t1+2*t3), np.sqrt(2)/3*(t1-t3)),\
						(np.sqrt(2./3)*t2*np.exp(2j*phi),np.sqrt(2)/3*(t1-t3), 1./3*(2*t1+t3))] )

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

def energies(ham_2d, args, NX = 40, NY = 40):
	"energies(ham_2d, args,NX=40, NY=40"
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

def enplot(args, NX = 40, NY = 40):
	pxgrid, pygrid, es = energies(hamiltonian_3, args, NX, NY)
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.set_xlim([-math.pi,math.pi])
	ax.set_ylim([-math.pi,math.pi])
#	ax.set_zlim([0,0.5])
#	ax.plot_wireframe(pxgrid, pygrid, es[1]-es[0], color = 'red')
#	ax.contour(pxgrid, pygrid, es[2]-es[1], 20)

	ax.plot_wireframe(pxgrid, pygrid, es[0], color = 'blue')
	ax.plot_wireframe(pxgrid, pygrid, es[1], color = 'red')
	ax.plot_wireframe(pxgrid, pygrid, es[2], color = 'green')

#	ax.contour(pxgrid, es[0], pygrid, 30)
#	ax.contour(pxgrid, es[1], pygrid, 30)
#	ax.contour(pxgrid, es[2], pygrid, 30)

#	fig.savefig('projection_levels.jpg')
	plt.show()

enplot((4,0.1,.05,.05,.1,.05,.05))
