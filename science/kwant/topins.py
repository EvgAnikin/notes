import kwant
import matplotlib.pyplot as plt
import numpy as np
import tinyarray

lat = kwant.lattice.general([[1,0], [0,1]])
symmetry = kwant.TranslationalSymmetry((0,1))
stripe = kwant.Builder(symmetry)

NX = 100
a = 5
A, B, D, M = 364.5, -686, 0, -10

xi = M
t = A/a
m1 = -0.25*a**2/(B + D)
m2 = -0.25*a**2/(B - D)

print(xi, 1/m1, 1/m2, t)

for i in range(NX):
    stripe[lat(i, 0)] = tinyarray.array([[xi + 2/m1, 0],
                                         [0, -xi - 2/m2]])
    stripe[lat(i, 0), lat(i, 1)] = tinyarray.array([[-0.5/m1, -t],
                                                    [t, 0.5/m2]])

for i in range(NX-1):
    stripe[lat(i, 0), lat(i+1, 0)] = tinyarray.array([[-0.5/m1,  -1j*t],
                                                      [-1j*t, 0.5/m2]])

stripe = stripe.finalized()
kwant.plotter.bands(stripe, momenta=np.linspace(-0.5, 0.5, 100))
