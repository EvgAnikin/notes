import kwant    
from matplotlib import pyplot

def stadium(position):
    x, y = position
    x = max (abs(x)-7, 0)
    return x**2 + y**2 < 10**2

sys = kwant.Builder()
sqlat = kwant.lattice.square()

sys[sqlat.shape(stadium, (0, 0))] = 4
sys[sqlat.neighbors()] = -1

lead_symmetry = kwant.TranslationalSymmetry([0, -1])
for start, end in [(-9, -3), (4, 10)]:
    lead = kwant.Builder(lead_symmetry)
    lead[(sqlat(x, 0) for x in range(start, end))] = 4
    lead[sqlat.neighbors()] = -1
    sys.attach_lead(lead)

sys = sys.finalized()

energies = [0.5 + 1e-4*i for i in range(300)]
conductances = [kwant.smatrix(sys, en).transmission(1, 0)
                                        for en in energies]
local_dos = kwant.ldos(sys, energy=0.4)
pyplot.plot(energies, conductances)
pyplot.show()
kwant.plotter.map(sys, local_dos, num_lead_cells=10)
