import math
import matplotlib.pyplot as plt
import numpy as np

input_file = open('gf_array', 'r')
np_file = np.load(input_file)

omegas = np_file['omegas']
gf = np_file['gf']

plt.plot(gf, omegas)
plt.plot(-gf[::-1], -omegas[::-1])
omega_range = np.array([min(omegas), max(omegas)])
#plt.plot(omega_range*2, np.ones(2)*0.05)
plt.ylim(tuple(omega_range*1.01))
plt.show()
