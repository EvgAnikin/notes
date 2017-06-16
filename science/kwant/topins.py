import kwant
import matplotlib.pyplot as plt
import tinyarray

lat = kwant.lattice.general([[1,0], [0,1]])
symmetry = kwant.TranslationalSymmetry((0,1))
stripe = kwant.Builder(symmetry)

NX = 20
xi, m1, m2, t = 0.1, 1, 1.5, 0.4 

for i in range(NX):
    stripe[lat(i, 0)] = tinyarray.array([[xi + 2/m1, 0],
                                         [0, -xi - 2/m2]])
    stripe[lat(i, 0), lat(i, 1)] = tinyarray.array([[-0.5/m1, -t],
                                                    [t, 0.5/m2]])

for i in range(NX-1):
    stripe[lat(i, 0), lat(i+1, 0)] = tinyarray.array([[-0.5/m1, -1j*t],
                                                         [-1j*t,  0.5/m2]])

stripe = stripe.finalized()
kwant.plotter.bands(stripe)
