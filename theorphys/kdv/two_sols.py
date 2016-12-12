import math
import numpy as np
import matplotlib.pyplot as plt

def two_sols(x, ka1, ka2, b1, b2):
	t11 = (ka1+ka2)/(ka1 - ka2)*b1
	t12 = 2*ka1/(ka1 - ka2)*b1
	t21 = 2*ka2/(ka2 - ka1)*b2
	t22 = (ka1+ka2)/(ka2 - ka1)*b2

	detA = 1 + t11*np.exp(-2*ka1*x) + t22*np.exp(-2*ka2*x) \
		+(t11*t22 - t12*t21)*np.exp(-2*(ka1+ka2)*x)
	ddetA = -2*ka1*t11*np.exp(-2*ka1*x) - 2*ka2*t22*np.exp(-2*ka2*x) \
		-2*(ka1+ka2)*(t11*t22-t12*t21)*np.exp(-2*(ka1+ka2)*x)
	d2detA = 4*ka1**2*t11*np.exp(-2*ka1*x) + 4*ka2**2*t22*np.exp(-2*ka2*x) \
		+4*(ka1+ka2)**2*(t11*t22-t12*t21)*np.exp(-2*(ka1+ka2)*x)

	return d2detA/detA - (ddetA/detA)**2

#x = np.linspace(-20,20,500)
#
#plt.plot(x, two_sols(x,0.1,0.2,1,3), 'r-')
#plt.show()
