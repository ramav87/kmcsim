#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import os
import re
import numpy as np
#import model
#import io
#from . import model
#from . import io
from .model import KMCModel
from .io import read_cfg, read_pars, write_cfg

class RunSim:

    # simulation flow control parameters with defaults
    t_max = 100.0
    print_period = 1.0
    save_traj_period = 10.
    measure_period = t_max
    param_file = 'params.in'
    incfg_file = 'init.xyz'
    outcfg_file = 'output.xyz'
    traj_file = 'kmc.trj'
    stats_file = 'statistics.dat'

    def __init__(self):
        pass

    def read(self, setup_file):
        """
        Read KMC input

        Parameters
        ----------
        setup_file : str
                     file containing information for setting up the simulation
        """

        path = os.path.dirname(setup_file)

        with open(setup_file, 'r') as f:

            # time limit (ms)
            self.t_max = float(re.findall('\S+', f.readline())[-1])

            # print runtime info period
            self.print_period = float(re.findall('\S+', f.readline())[-1])

            # save trajectory period
            self.save_traj_period = float(re.findall('\S+', f.readline())[-1])

            # kmx parameter file - rates (1/ms)
            self.param_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # input configuration file
            self.incfg_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # input configuration file
            self.outcfg_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # trajectory file
            self.traj_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # statistics file
            self.stats_file = os.path.join(path, re.findall('\S+', f.readline())[-1])


    def init_sim(self):
        """
        Initialize KMC simulation: build model and set its parameters
        """

        # read input configuration
        lat_type, box, xyz = read_cfg(self.incfg_file)

        # Initialize KMC system with appropriate lattice type
        self.kmc = KMCModel(lat_type)
        self.kmc.make_lattice(xyz, box)

        # read kMC parameters (reaction rates)
        rates = read_pars(self.param_file)

        # make event list (e.g., identify deposition sites)
        self.kmc.init_events(rates)



    def run(self, random_seed=42):
        """
        Run KMC simulation
        """

        # initial values
        t = t_print = t_save = t_measure = 0.0
        it = 0

        # print initial numbers
        if (self.print_period < self.t_max):
            print('time, iteration, number of atoms')
            print(t, it, self.kmc.nat)

        while t < self.t_max:

            self.kmc.step()

            t += self.kmc.advance_time()

            it += 1

            # perform runtime outputs
            if (t - t_print) > self.print_period:
                print(t, it, self.kmc.nat)
                t_print = t

            if (t - t_save) > self.save_traj_period:
                t_save = t

            if (t - t_measure) > self.measure_period:
                t_measure = t

        if (self.print_period < self.t_max):
            print('End of simulation')
            print(t, it, self.kmc.nat)


    def output(self):
        """
        Save the final state and statistics
        """

        xyz, box, grain = self.kmc.get_conf()
        write_cfg(self.outcfg_file, xyz, box, grain)


if __name__ == "__main__":

    sim = RunSim()

    sim.read(sys.argv[0])

    sim.init_sim()

    sim.run()

    sim.output()

