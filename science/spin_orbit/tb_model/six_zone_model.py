from __future__ import division
import math
import numpy as np

def valence_ham(tpar, tperp, px, py, pz):
    ham = np.diag([(tpar + tperp)*(math.cos(px) + math.cos(py)) + 2*tperp*math.cos(pz),
                   (tpar + 5*tperp)/3*(math.cos(px) + math.cos(py)) + 
