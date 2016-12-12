from en_levels import *

args = 0.14,0.2,0.3,0.15,0.15,1,0

e1 = linalg.eigvalsh(hamiltonian(*args))
e2 = linalg.eigvalsh(hamiltonian_2(*args))

print e1
print e2
