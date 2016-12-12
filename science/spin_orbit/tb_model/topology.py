import numpy as np
import matplotlib.pyplot as plt

def T(phi, *args): 
	t1,t2,t3 = args
	return  np.array( [(t1, 1/np.sqrt(3)*t2*np.exp(-2j*phi),\
									 np.sqrt(2./3)*t2*np.exp(-2j*phi)),\
						(1/np.sqrt(3)*t2*np.exp(2j*phi), 1./3*(t1+2*t3), np.sqrt(2)/3*(t1-t3)),\
						(np.sqrt(2./3)*t2*np.exp(2j*phi),np.sqrt(2)/3*(t1-t3), 1./3*(2*t1+t3))] )

def ham(px,py, ESO,t1,t2,t3,td1,td2,td3):
	"hamiltonian_3(px,py,ESO,t1,t2,t3,td1,td2,td3)"
	atomic = np.diag([0,0,-ESO])	

	return atomic + 2*np.cos(px)*T(0,t1,t2,t3) + 2*np.cos(py)*T(np.pi/2,t1,t2,t3) +\
					2*np.cos(px+py)*T(np.pi/4,td1,td2,td3) + 2*np.cos(px-py)*T(-np.pi/4,td1,td2,td3)

def a(px,py, *args):
	return ham(px,py,*args)[0,0]

def b(px,py, *args):
	return ham(px,py,*args)[0,1]

def c(px,py, *args):
	return ham(px,py,*args)[1,1]

def e(px,py, *args):
	A = a(px,py,*args)
	B = b(px,py,*args)
	C = c(px,py,*args)
	return 0.5*(A + C - np.sqrt((A-C)**2 + 4*abs(B)**2))

def phase_mult(px,py,*args):
	A = a(px,py,*args)
	B = b(px,py,*args)
	C = c(px,py,*args)
	E = e(px,py,*args)

	return (E - C)/B * np.sqrt( ((E - A)**2 + abs(B)**2)/((E - C)**2 + abs(B)**2) )

p = 0.01
args = 2,0.1,0.05,0.05,0.1,0.05,0.05

phis = np.linspace(0,2*np.pi,100)
mults_re = []
mults_im = []

for phi in phis:
	mults_re.append(phase_mult(p*np.cos(phi), p*np.sin(phi), *args).real)
	mults_im.append(phase_mult(p*np.cos(phi), p*np.sin(phi), *args).imag)


plt.plot(phis,mults_re)
plt.plot(phis,mults_im)
plt.xlim(0,2*np.pi)
plt.ylim(-1.2,1.2)
plt.show()
