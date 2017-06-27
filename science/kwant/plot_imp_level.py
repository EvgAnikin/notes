import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

input_file = open('gf_array', 'r')
np_file = np.load(input_file)

omegas = np_file['omegas']
gf = np_file['gf']

interpolated_omega = interp1d(1/gf, omegas, kind='cubic')
x = np.linspace(np.amin(1/gf), np.amax(1/gf), 50)
plt.plot(x, interpolated_omega(x))

#plt.plot(gf, omegas)
#plt.plot(-gf[::-1], -omegas[::-1])
omega_range = np.array([min(omegas), max(omegas)])
#plt.plot(omega_range*2, np.ones(2)*0.05)
plt.ylim(tuple(omega_range*1.01))
plt.show()
