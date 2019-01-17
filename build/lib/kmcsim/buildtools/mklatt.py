#!//anaconda/envs/py36/bin/python
#
# File name:   mklatt.py
# Date:        2018/08/02 16:15
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re

def make_sc(box):
    print('sc', *box)
    return

def make_bcc(box):
    print('bcc', *box)
    return

def make_fcc(box):

    lx, ly, lz = box

    latt = {}
    latt['nat'] = lx*ly*lz
    latt['box'] = ['fcc', 2*lx, ly, lz]
    latt['xyzs'] = []

    # box dimensions in lattice units

    # layer number
    for iz in range(lz):
        # layer structure
        for iy in range(ly):
            for ix in range(lx):
                rx = 2*ix + (iy + iz)%2
                latt['xyzs'].append(['Ni', rx, iy, iz])

    return latt

def write_latt(latt, fname):

    with open(fname, 'w') as fo:
        fo.write('{0:d}\n'.format(latt['nat']))
        fo.write('{0} {1:d} {2:d} {3:d}\n'.format(*latt['box']))
        for row in latt['xyzs']:
            fo.write('{0} {1:d} {2:d} {3:d}\n'.format(*row))


if __name__ == "__main__":

    # dictionary of lattice build functions
    make_lattice = {'fcc':make_fcc, 'bcc':make_bcc, 'sc':make_sc}

    # read lattice type and parameters
    with open(sys.argv[1], 'r') as f:

        # lattice type
        ltype = re.findall('\S+', f.readline())[0]

        # box dimensions
        box = [int(d) for d in re.findall('\S+', f.readline())]

    # make lattice
    latt = make_lattice[ltype](box)

    write_latt(latt, 'init.xyz')

# end of mklatt.py 
