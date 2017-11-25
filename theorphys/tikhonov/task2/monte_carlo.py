from __future__ import division
import math
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad 
from scipy.optimize import fsolve 


def metropolis_step(X, alpha, tau):
    N = len(X)
    i = np.random.randint(N)

    def partial_action(x):
        return (x - X[i-1])**2/2/tau + (x - X[(i+1)%N])**2/2/tau + (x**2/2 + alpha/24*x**4)*tau

    x_max = 8
    x = -x_max + 2*x_max*np.random.rand()

    old_weight = math.exp(-partial_action(X[i]))
    new_weight = math.exp(-partial_action(x))

    X_new = X.copy()
    X_new[i] = x

    p = np.random.rand()
    return X_new if p < new_weight/old_weight else X


N = 15
beta = 2
tau = beta/N
alpha = 0

#X = np.random.rand(N) - 0.5*np.ones(N)
X = np.zeros(N)
x_squares = []

N_samples = 100
N_interval = 100000
for i in xrange(N_samples):
    for j in xrange(N_interval):
        X = metropolis_step(X, alpha, tau)
    x_squares.append(np.mean(X**2))

print np.mean(x_squares), np.std(x_squares)/math.sqrt(N_samples)
