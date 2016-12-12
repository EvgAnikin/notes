import math
import matplotlib.pyplot as plt
import numpy as np

input_file = open('gf_array', 'r')
np_file = np.load(input_file)

omegas = np_file['arr_0']
gf = np_file['arr_1']

plt.plot(omegas, gf)
plt.plot(omegas, -gf[::-1])
plt.show()
