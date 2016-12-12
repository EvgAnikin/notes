import sys
import math
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

def hamiltonian(n, py, *args):
	"hamiltonian(n,py,*args), args : ESO, t1,t2,t3"
	ESO, t1,t2,t3 = args
	tlong = t1 + t2
	ttrans = t1 - t2

	atom_ham =-ESO/3. * np.array([(1,1j,-1), \
							(-1j,1,1j), \
							(-1, -1j, 1.)])

	stripe = linalg.block_diag(*(n*[atom_ham]))

	for i in xrange(n):
		stripe[3*i,3*i]    += 2*ttrans*math.cos(py)
		stripe[3*i+1,3*i+1]+= 2*tlong *math.cos(py)
		stripe[3*i+2,3*i+2]+= 2*t3*math.cos(py)

	for i in xrange(n-1):
		stripe[3*i, 3*(i+1)] = tlong
		stripe[3*(i+1), 3*i] = tlong

		stripe[3*i+1, 3*(i+1)+1] = ttrans
		stripe[3*(i+1)+1, 3*i+1] = ttrans

		stripe[3*i+2, 3*(i+1)+2] = t3
		stripe[3*(i+1)+2, 3*i+2] = t3
	
	return stripe

def T(phi, *args): 
	"T(phi,t1,t2,t3)"
	t1,t2,t3 = args
	return  np.array( [(t1, 1/np.sqrt(3)*t2*np.exp(-2j*phi),\
									 np.sqrt(2./3)*t2*np.exp(-2j*phi)),\
						(1/np.sqrt(3)*t2*np.exp(2j*phi), 1./3*(t1+2*t3), np.sqrt(2)/3*(t1-t3)),\
						(np.sqrt(2./3)*t2*np.exp(2j*phi),np.sqrt(2)/3*(t1-t3), 1./3*(2*t1+t3))] )

def hamiltonian_1(n,py,*args):
	"hamiltonian(n,py,*args), args : ESO, t1,t2,t3,td1,td2,td3"
	ESO,t1,t2,t3,td1,td2,td3 = args
	atom_ham = np.diag([0,0,-ESO])

	stripe = linalg.block_diag(*([atom_ham + 2*np.cos(py)*T(np.pi/2,t1,t2,t3)]*n))
	
	for i in xrange(0,n-1):
		stripe[3*i:3*(i+1),3*(i+1):3*(i+2)] =\
				 T(0,t1,t2,t3)\
				+np.exp(-1j*py)*T( np.pi/4,td1,td2,td3)\
				+np.exp( 1j*py)*T(-np.pi/4,td1,td2,td3)

	 	stripe[3*(i+1):3*(i+2),3*i:3*(i+1)] =\
				T(0,t1,t2,t3)\
				+np.exp( 1j*py)*T( np.pi/4,td1,td2,td3)\
				+np.exp(-1j*py)*T(-np.pi/4,td1,td2,td3)

	return stripe

	
	

if __name__ == '__main__':

	if len(sys.argv[1:]) == 7:
		ESO,t1,t2,t3,td1,td2,td3 = map(float, sys.argv[1:])	
	else:
		print 'Usage : ESO,t1,t2,t3,td1,td2,td3'
		exit()

	args = ESO,t1,t2,t3,td1,td2,td3

	N = 60
	NP = 200
	
#	allpy = np.linspace(1.8,2.4, NP)
	allpy = np.linspace(0,2*np.pi, NP)
	levels = []
	for py in allpy:
		levels.append(linalg.eigvalsh(hamiltonian_1(N,py,*args)))

	for energies in zip(*levels):
		plt.plot(allpy, energies, 'r-')
	plt.xlim((0,2*np.pi))
#	plt.ylim((-1.5,0.5))
#	plt.savefig('stripe.jpg')
	plt.show()
	
