from __future__ import division
import math
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad 
from scipy.optimize import fsolve 
import sys


def metropolis_step(X, alpha, tau):
    N = len(X)
    i = np.random.randint(N)

    def partial_action(x):
        return (x - X[i-1])**2/2/tau + (x - X[(i+1)%N])**2/2/tau#+ (x**2/2 + alpha/24*x**4)*tau

    x_max = 8
    x = -x_max + 2*x_max*np.random.rand()

    old_weight = math.exp(-partial_action(X[i]))
    new_weight = math.exp(-partial_action(x))

    p = np.random.rand()
    if p < new_weight/old_weight:
        X[i] = x


def heatbath_step(X, tau):
    N = len(X)
    i = np.random.randint(N)

    X[i] = np.random.normal(loc=(X[i-1] + X[(i+1)%N])/(2 + tau**2), 
                            scale=math.sqrt(tau/(2 + tau**2)))


N_samples = int(sys.argv[1])
N_interval = int(sys.argv[2])
N = int(sys.argv[3])
beta = float(sys.argv[4])
tau = beta/N
alpha = float(sys.argv[5])

#X = np.random.rand(N) - 0.5*np.ones(N)
X = np.zeros(N)
x_squares = []
correlation = []

N_initial = 20000

for i in xrange(N_initial):
    heatbath_step(X, tau)

Xarray = np.zeros((N, N_samples*N_interval))

for i in xrange(N_samples):
    for j in xrange(N_interval):
        heatbath_step(X, tau)
    Xarray[:,i] = X


correlation = np.zeros(N_samples)
for i in xrange(N_samples):
    correlation[i] = np.mean([np.sum(Xarray[:,j]*Xarray[:,j+i]) 
                              for j in xrange(0, N_samples - j)])


plt.plot(range(0, N_samples*N_interval, N_interval), correlation)
plt.show()

#print np.mean(x_squares), np.std(x_squares)/math.sqrt(N_samples)
