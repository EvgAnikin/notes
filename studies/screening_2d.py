import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import *
from scipy.special import *

def integrand(x, x0):
    return jn(0, x)/(x + x0)
    

def screened_coulomb(x0):
    return 1/x0 - quad(integrand, 0, 1000, args = (x0), limit=1000)[0]


x = np.linspace(0.1,100,200)

potential = []
for xi in x:
    potential.append(screened_coulomb(xi))

plt.plot(x, np.log(potential))
plt.show()
