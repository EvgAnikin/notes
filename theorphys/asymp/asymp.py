import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
 
def fact(n):
	res = 1
	for i in range(1,n+1):
		res *= i
	return res

def approx_Efunc(x,n):
	res = 0
	for k in range(0,n+1):
		res += (-1)**k * fact(k) / x**(k+1)
	return res

def series_Efunc(x,n):
	res = -np.log(x) + integrate.quad(lambda t : (np.exp(-t) - 1)/t, 0, 1)[0] + \
							integrate.quad(lambda t : np.exp(-t)/t, 1, np.inf)[0]
	for k in range(1,n+1):
		res += (-1)**(k+1)*x**k/k/fact(k)
	return res*np.exp(x)

def Efunc(x):
	return integrate.quad(lambda t : np.exp(-x*t)/(1 + t), 0, np.inf)[0]

NPOINTS = 100
XMAX = 2
NAPPROX = 5
NSERIES = 3
epsilon = 0.001

fig = plt.figure()
ax = fig.add_subplot(111)

x = np.linspace(0+epsilon,XMAX, NPOINTS)

efunc_val = []
for xi in x:
	efunc_val.append(Efunc(xi))

efunc_val = np.array(efunc_val)

ax.plot(x, efunc_val, color = 'black', linewidth = 2.0, label = 'function')
#ax.set_yscale("log")
plt.xlim(0, XMAX)
#plt.ylim(1e-6, 1e1)
plt.ylim(0,10)
#for n in xrange(1,NAPPROX+1):
#	ax.plot(x, abs((approx_Efunc(x,n) - efunc_val)/efunc_val), label = '%d th approximation' % n)
#
for n in xrange(0,NSERIES+1):
	ax.plot(x, series_Efunc(x,n), label = '%d th series' % n)

ax.legend()
plt.show()
