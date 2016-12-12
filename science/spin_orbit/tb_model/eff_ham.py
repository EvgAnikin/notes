from tblib import *
import matplotlib.pyplot as plt
from en_levels import *

def effective_hamiltonian(p,ESO, m1,m2):
    H0 = np.diag( [p**2/m1/2, p**2/m2/2, p**2/m2/2] )
    
    V =	-ESO/3. * np.array([(1,1j,-1), \
						(-1j,1,1j), \
						(-1, -1j, 1.)])
	
    return H0 + V

if __name__ == '__main__':
#    args = 0.2,-1,-2
    args = 0.2,0.15,0.5,0.3,0.15,0.1,1,0
    h_mod = lambda *args : s_type_hamiltonian(0, *args)
    p,E = stripe_energies(h_mod, args, plim = (-1,1))
    
    for level in E:
        plt.plot(p,level)

    plt.xlim(-1,1)
#    plt.ylim(-1,0.1)
    plt.show()
