#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek, Rama Vasudevan
#
# Description: 
#

import sys
import os
import re
import numpy as np
import time
import random
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
    save_traj_period = 10.0
    measure_period = t_max
    param_file = 'params.in'
    incfg_file = 'init.xyz'
    outcfg_file = 'output.xyz'
    traj_file = 'kmc.trj'
    stats_file = 'statistics.dat'

    def __init__(self, random_seed=42):
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
            self.t_max = float(re.findall('\S+', f.readline())[-1]) #time for whole simulation

            self.max_time_steps = float(re.findall('\S+', f.readline())[-1]) #maximum time steps

            # print runtime info period
            self.print_period = float(re.findall('\S+', f.readline())[-1])

            # print runtime info period
            self.update_rate_period = float(re.findall('\S+', f.readline())[-1]) #time to update rate

            # save trajectory period
            self.save_traj_period = float(re.findall('\S+', f.readline())[-1])

            # kmc parameter file - rates (1/ms)
            self.param_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # input configuration file
            self.incfg_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # output configuration file
            self.outcfg_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # trajectory file
            self.traj_file = os.path.join(path, re.findall('\S+', f.readline())[-1])

            # statistics file
            self.stats_file = os.path.join(path, re.findall('\S+', f.readline())[-1])


    def init_sim(self, random_seed=42):
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
        
        # initialize random number generator
        random.seed(random_seed)
        
        #setup global time step
        self.global_time_steps =0
        self.global_time = 0

        #setup an output file counter
        self.output_count = 0

    def run(self):
        """
        Run KMC simulation
        Verbose: False (Default). Set to True for event prints
        """
        self.kmc.verbose = False
        # initial values
        t = t_print = t_save = t_measure = 0.0
        it = 0

        # print initial numbers
        if (self.print_period < self.t_max):
            print('time, iteration, number of atoms')
            print(t, it, self.kmc.nat)

        while t < self.t_max: #e.g. while t<max_iterations (it_max..).
        #while it<self.max_time_steps:
            t += self.kmc.advance_time()
            self.global_time = t
            self.kmc.step()

            it += 1

            # perform runtime outputs and any rate changes
            if (t - t_print) > self.print_period:
                print(t, it, self.kmc.nat)
                t_print = t
                xyz, box, grain, atom_types = self.kmc.get_conf()
                self.output_count += 1
                if self.output_count<10:
                    output_filename = self.outcfg_file[:-4] + '0' + str(self.output_count) + '.xyz'
                write_cfg(output_filename, xyz, box, grain, atom_types)
           
            if (t - t_save) > self.save_traj_period:
                t_save = t

                print("save trajectory", t, it, self.kmc.nat)

            if (t - t_measure) > self.measure_period:
                t_measure = t

        if (self.print_period < self.t_max):
            print('End of simulation')
            print(t, it, self.kmc.nat)


    def run_to_next_step(self, random_seed = 0.42):
        #run the simulation only to the next step
        if self.global_time>= self.t_max:
            #print("Time steps {} and Max time steps {}".format(self.global_time_steps, self.max_time_steps))
            #print('End of simulation')
            return 1
        else:
            it=0
            t_local=0
            global_time_start = self.global_time
            while t_local<self.update_rate_period:
                t_local += self.kmc.advance_time()
                self.kmc.step()
                it += 1
                self.global_time_steps += 1 #update the global counter.
            self.global_time =global_time_start + t_local
            return 0
            #print("iterations {}".format(it))
            
    def output(self):
        """
        Save the final state and statistics
        """

        xyz, box, grain, atom_types = self.kmc.get_conf()
        self.output_count+=1
        output_filename = self.outcfg_file[:-4] + str(self.output_count) + '.xyz'
        write_cfg(output_filename, xyz, box, grain,atom_types )

    def update_rate(self, new_rates, verbose=True):
        if verbose:
            print('current rates:{} '.format(self.kmc.etree.rates ))
            print('new rates:{}'.format(new_rates))
            print('kmc step:', self.global_time_steps)
            
        #Update global rates
        self.kmc.etree.update_global_rate(new_rates)
        #Update events
        n_events = np.array([len(e) for e in self.kmc.event_list])
        self.kmc.etree.update_events(n_events)

if __name__ == "__main__":

    sim = RunSim()

    sim.read(sys.argv[0])

    sim.init_sim()

    sim.run()

    sim.output()

