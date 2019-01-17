#
# File name:   runsim.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
from .runsim import RunSim

if __name__ == "__main__":

    sim.RunSim()

    sim.read(sys.argv[0])

    sim.initialize()

    sim.run()

    sim.output()

