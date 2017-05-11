import numpy as np
from lattice_2d import draw_state
import sys

densities = np.load(sys.argv[1])['densities']
