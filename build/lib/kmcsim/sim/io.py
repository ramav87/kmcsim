#!//anaconda/envs/py36/bin/python
#
# File name:   io.py
# Date:        2018/08/16 21:56
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np

def read_cfg(file_name):
    """Read input configuration and create simulation box"""

    with open(file_name, 'r') as f:

        # number of atoms
        nat = int(re.findall('\S+', f.readline())[0])

        sarr = re.findall('\S+', f.readline())

        # lattice type
        lat_type = sarr[0]

        # simulation box size
        box = np.array([float(l) for l in sarr[1:]], dtype=int)

        # read particle coordinates and grain identity
        xyz = []
        for i in range(nat):
            xyz.append(np.array([int(l) for l in re.findall('\S+', f.readline())[1:]]))

    return lat_type, box, xyz

def write_cfg(file_name, xyz, box, grain):
    """Write output configuration to xyz file"""

    with open(file_name, 'w') as f:
        # number of atoms and box
        nat = len(xyz)
        f.write(f'{nat}\n')
        f.write(f'fcc {box[0]} {box[1]} {box[2]}\n')

        # write particle coordinates and grain identity
        for r, g in zip(xyz, grain):
            f.write(f'{g} {r[0]} {r[1]} {r[2]}\n')

def read_pars(file_name):
    """Read rates for different event types"""

    with open(file_name, 'r') as f:
        rates = []
        for line in iter(f.readline, ''):
            rates.append(float(re.findall('\S+', line)[1]))

    return np.array(rates)

def read_input(file_name):
    """
    First reads top level simulation input file
    and then configuration and parameter files.

    Parameters
    ----------
    file_name : str
                input file name

    Returns
    -------
    input_data : dict
                 dictionary with configuration and rates data
    """

    with open(file_name, 'r') as f:
        # time limit (ms)
        t_max = float(re.findall('\S+', f.readline())[1])

        # input configuration file
        incfg_file = re.findall('\S+', f.readline())[1]

        # kmx parameter file - rates (1/ms)
        pars_file = re.findall('\S+', f.readline())[1]

    # read configuratio and parameters (rates) file
    xyz, box = read_cfg(incfg_file)
    rates = read_pars(pars_file)

    # construct dictionary of the data
    input_data = {'t_max':t_max}
    input_data['xyz'] = xyz
    input_data['box'] = box
    input_data['rates'] = rates

    return input_data

