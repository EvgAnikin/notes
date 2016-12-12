import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def solution(init_func, nsteps, tau, a):
	result = []

	prev_slice = init_func
	for i in xrange(nsteps):
		n = int(2.*tau/a**2) + 1	
		for j in xrange(n):
			next_slice = prev_slice + tau/n/a**2*(np.roll(prev_slice,1) +\
				np.roll(prev_slice,-1) - 2*prev_slice)
			prev_slice = next_slice
		result.append(next_slice)
	return np.array(result)

N = 100
NSTEPS = 100
XMAX = 25
TMAX = 3
x0 = 0.2

x = np.linspace(-XMAX,XMAX,N)
init = 1./2./np.cosh(x/2./x0)

a = 2.*XMAX/N
tau = float(TMAX)/NSTEPS
t = [tau*i for i in xrange(1,NSTEPS+1)]

sol = solution(init,NSTEPS,tau, a)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

xgrid,tgrid = np.meshgrid(x,t)
start = 20
ax.plot_surface(xgrid[start:],tgrid[start:], \
			(sol*np.cos(np.pi*xgrid/8/tgrid)/np.exp(-xgrid**2/4/tgrid**2))[start:])
ax.set_zlim([0,2])
plt.show()
#print np.exp(-xgrid**2/4./tgrid**2)
