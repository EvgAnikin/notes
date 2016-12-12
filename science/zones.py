# coding:utf-8

# Программа строит зависимость энергии частицы в поле потенциальной 
# гребёнки от квазиимпульса.

from math import *
import matplotlib.pyplot as plt

kappa = -3

def qp(p0):
	arg = cosh(p0) + kappa/p0*sinh(p0)
	if -1 <= arg < 1:
		return acos(arg)
	else:
		return None

dp = 0.01
p0 = dp
pmax = 10
n = 0

quasi_mom = []
energy = []

while p0 < pmax:
	while qp(p0) == None and p0 < pmax:
		p0 += dp
	while qp(p0) and p0 < pmax:
		energy.append(-p0**2/2)
		energy.append(-p0**2/2)
		quasi_mom.append(   (n+1)/2 * (2*pi) + (-1)**n * qp(p0))
		quasi_mom.append( -((n+1)/2 * (2*pi) + (-1)**n * qp(p0)))
		p0 += dp
	n += 1

for tup in zip(quasi_mom, energy):
	print tup[0], tup[1]
