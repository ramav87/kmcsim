import sys
import os
import numpy as np

from kmcsim.buildtools import make_fcc, write_latt
from kmcsim.sim import KMCModel
from kmcsim.sim import EventTree
from kmcsim.sim import RunSim

wdir = r'data/working'

# make substrate: perfect FCC lattice of given dimensions
box = [16, 16, 4]
latt = make_fcc(box)

# extend the box in the z-direction to make space for new layers to grow
latt['box'][3] = 256

# write initial configuration to xyz file in the working directory
write_latt(latt, os.path.join(wdir,'ni.xyz'))

sim = RunSim()

sim.read(os.path.join(wdir,'kmc.input'))

sim.init_sim()

sim.run()

sim.output()
