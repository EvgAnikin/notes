import numpy as np
import matplotlib.pyplot as plt

def next_psi(psi, psi_prev, e, delta_e):
    psi_next = (2 - e + delta_e)*psi - psi_prev
    return psi_next, psi


def iterate(e, magnitude, n):
    wf = []
    psi, psi_prev = 1, 0
    for i in xrange(n):
        delta_e = (-1 + 2*np.random.rand())*magnitude
        wf.append(psi)
        psi, psi_prev = next_psi(psi, psi_prev, e, delta_e)
    wf.append(psi)
    return wf


if __name__ == '__main__':
    wf = np.array(iterate(0.2, 0.02, 1000000))
    plt.plot(np.log(abs(wf[::1000])))
    plt.show()
