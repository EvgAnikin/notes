import math
import matplotlib.pyplot as plt
import numpy as np

def angle_potential(theta, alpha=math.pi/2):
    return 1/(np.sin(math.pi*theta/(2*math.pi - alpha)))**(2 - alpha/math.pi)

ax = plt.subplot(111, projection='polar')

alpha = math.pi/180.*45
nlines = 40
epsilon = 0.02
theta = np.linspace(epsilon, 2*math.pi - alpha - epsilon, 400)
delta_r = 0.1
rmax = 10

for i in xrange(1, nlines+1):
    ax.plot(theta, ((i*delta_r)**(2 - alpha/math.pi)*
                    angle_potential(theta, alpha)), color='red')

ax.plot([0,0], [0,rmax], color='black')
ax.plot([2*math.pi - alpha, 2*math.pi - alpha], [0, rmax], color='black')
ax.grid(False)

ax.set_rmax(rmax)
plt.show()
