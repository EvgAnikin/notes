from __future__ import division
import math
import numpy as np
from tblib import stripe_energies
import matplotlib.pyplot as plt

def valence_ham(px, py, pz, tpar, tperp, Ev):
    ham =  ((Ev - (2*tpar + 4*tperp))*np.identity(4) + 
           np.diag([(tpar + tperp)*(math.cos(px) + math.cos(py)) + 2*tperp*math.cos(pz),
                   ((tpar + 5*tperp)/3*(math.cos(px) + math.cos(py)) +
                            (2*tperp/3 + 4*tpar/3)*math.cos(pz)),
                   ((tpar + 5*tperp)/3*(math.cos(px) + math.cos(py)) +
                            (2*tperp/3 + 4*tpar/3)*math.cos(pz)),
                   (tpar + tperp)*(math.cos(px) + math.cos(py)) + 2*tperp*math.cos(pz)]))
    ham[0,2] = ham[2,0] = ham[1,3] = ham[3,1] = (-1/math.sqrt(3)*(tpar - tperp)*
                                                        (math.cos(px) - math.cos(py)))
    return ham

def kane_tb_ham(px, py, pz, args):
    ts, tpar, tperp, P, Ec, Ev = args
    ham = np.zeros(6*6, dtype=complex).reshape(6,6)
    ham[2:,2:] = valence_ham(px, py, pz, tpar, tperp, Ev)

    Ax = np.array([[-1/math.sqrt(2), 0, 1/math.sqrt(6), 0],
                   [0, -1/math.sqrt(6), 0, 1/math.sqrt(2)]])
    Ay = np.array([[-1j/math.sqrt(2), 0, -1j/math.sqrt(6), 0],
                   [0, -1j/math.sqrt(6), 0, -1j/math.sqrt(2)]])
    Az = np.array([[0, math.sqrt(2/3), 0, 0],
                   [0, 0, math.sqrt(2/3), 0]])
    ham[0:2,2:] = P*(math.sin(px)*Ax + math.sin(py)*Ay + math.sin(pz)*Az)
    ham[2:,0:2] = ham[0:2,2:].T

    ham[0:2,0:2] = np.identity(2)*(Ec + 2*ts*(3 - math.cos(px) - math.cos(py) - math.cos(pz)))

    return ham


hbar = 1.0545718e-27
me = 9.10938356e-28
eV = 1.6021766209e-12
a = 6.46e-8
E0 = hbar**2/2./me/a**2/eV

#gamma1 = 1.47
#gamma2 = -0.28
#Ep = 18.8
#Ec = 1.606 - 0.57
#Ev = -0.57

gamma1 = 4.1
gamma2 = 0.5
Ep = 18.8
Ec = -0.303
Ev = 0

tpar = E0*(gamma1 + 4*gamma2)
tperp = E0*(gamma1 - 2*gamma2)
ts = E0
P = math.sqrt(Ep*E0)

ham_x = lambda px: kane_tb_ham(px, 0., 0, (ts, tpar, tperp, P, Ec, Ev))

pyrange, energies = stripe_energies(ham_x, (), NX = 101, plim = (-math.pi,math.pi))

for band in energies:
    plt.plot(pyrange, band)

plt.xlim(-math.pi, math.pi)
plt.ylim(-1.5,1.5)
plt.show()
