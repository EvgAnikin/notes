import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation
from two_sols import *

XMIN = -40
XMAX = 60
DT = 0.1
kappa1 = 1.1
kappa2 = 1.0
b1 = 6e-23
b2 = -1e-16

fig = plt.figure()
ax = plt.axes(xlim=(XMIN, XMAX), ylim=(-2,2))
line, = ax.plot([], [], lw=2)

def init():
	line.set_data([],[])
	return line,

def animate(i):
	x = np.linspace(XMIN,XMAX,400)
	y = two_sols(x, kappa1,kappa2, b1*math.exp(8*kappa1**3*DT*i), b2*math.exp(8*kappa2**3*DT*i))
	print b1*math.exp(8*kappa1**3*DT*i), b2*math.exp(8*kappa2**3*DT*i)
	line.set_data(x,y)
	return line,

anim = animation.FuncAnimation(fig, animate, 
		init_func=init, frames=200, interval=25, blit=True, repeat=False)
anim.save('two_sols_animation.mp4', writer='mencoder', fps=25)
