import math
import matplotlib.pyplot as plt
import numpy as np

def deltaE(omega, dt, ksi, t):
	chk = abs((omega**2 - ksi**2)/(2*t**2) - 1)
	if abs(omega) >= math.sqrt(ksi**2 + 4*t**2):
		return omega - ksi - 2*(1 + chk - math.sqrt(chk**2 - 1))/(omega + ksi)*(t + dt)**2
	elif abs(omega) <= abs(ksi):
		return omega - ksi - 2*(1 - chk + math.sqrt(chk**2 - 1))/(omega + ksi)*(t + dt)**2
	else:
		return None

ksi = 1
t   = 1.5
edge = math.sqrt(ksi**2 + 4*t**2)
dtmin = -5
dtmax = 5
N = 100

def plotfamily(omegas, color):
	for omega in omegas:
		args = np.linspace(dtmin, dtmax, N)
		plt.plot(args, map(lambda x : deltaE(omega, x, ksi,t), args), color + '-')

plotfamily(np.append((ksi*np.sin(\
		np.linspace(np.pi/2,-np.pi/2, 15, endpoint = False))), -ksi+0.001), 'r')
plotfamily(np.linspace(edge, 5*edge, 20), 'g')
plotfamily(np.linspace(-5*edge, -edge, 20), 'b')
plt.axis([dtmin, dtmax, -10, 10])
plt.savefig('levels_family.png')
