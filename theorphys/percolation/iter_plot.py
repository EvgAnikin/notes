#encoding:utf-8
#русский комментарий
import math
import matplotlib.pyplot as plt
import numpy as np


def f(x):
    return 2*x**2 - x**4


def plot_iteration(n, x0, ax):
    x = x0
    y = f(x)
    color = 'black'
    x_points = []
    y_points = []
    for i in xrange(n):
        x_points.extend([x, y])
        y_points.extend([y, y])
        x = y
        y = f(x)
    ax.plot(x_points, y_points, '-', color=color, label='iteration')

x = np.linspace(0,1,100)
y = f(x)


fig, ax = plt.subplots()
ax.plot(x,y, '-', color='b', label='$y = x^2 - x^4$')
ax.plot(x,x, '--', color='b', label='$y = x$')
x0 = (math.sqrt(5) - 1)/2 - 0.1
plot_iteration(10, x0, ax)

start_circle = plt.Circle((x0, f(x0)), 0.01, color='b', alpha=1, zorder=10)   
ax.add_artist(start_circle)
ax.legend(loc=2)

plt.xlim(0,1)
plt.ylim(0,1)

plt.show()
