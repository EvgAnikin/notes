import kwant
import matplotlib.pyplot as plt
import tinyarray

lat = kwant.lattice.general([[1,0], [0,1]])
symmetry = kwant.TranslationalSymmetry((0,1))
stripe = kwant.Builder(symmetry)

NX = 20
xi, m, t = -0.3, 1, 0.4 

for i in range(NX):
    stripe[lat(i, 0)] = tinyarray.array([[xi + 2/m, 0],
                                         [0, -xi - 2/m]])
    stripe[lat(i, 0), lat(i, 1)] = tinyarray.array([[-0.5/m, -t],
                                                    [t, 0.5/m]])

for i in range(NX-1):
    stripe[lat(i, 0), lat(i+1, 0)] = tinyarray.array([[-0.5/m, -1j*t],
                                                         [-1j*t,  0.5/m]])

stripe = stripe.finalized()
kwant.plotter.bands(stripe)
