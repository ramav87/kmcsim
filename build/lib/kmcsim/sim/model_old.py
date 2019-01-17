#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np
from itertools import product
from collections import Counter
#import events
from .events import EventTree

class KMCModel:
    """Class managing kmc moves and event modifications"""

    Rs = 1.0 # top node / sum of all rates

    def __init__(self, latt_type):
        self.latt_type = latt_type
        self.__setup_neighbors()

    def __setup_neighbors(self):
        """Create lists of neighboring (nn & nnn) sites"""

        nbrlist = []

        if self.latt_type == 'fcc':
            # NN
            nbrlist.append(np.array([ 1, 1, 0]))
            nbrlist.append(np.array([ 1,-1, 0]))
            nbrlist.append(np.array([-1, 1, 0]))
            nbrlist.append(np.array([-1,-1, 0]))
            nbrlist.append(np.array([ 1, 0,-1]))
            nbrlist.append(np.array([-1, 0,-1]))
            nbrlist.append(np.array([ 0, 1,-1]))
            nbrlist.append(np.array([ 0,-1,-1]))
            nbrlist.append(np.array([ 1, 0, 1]))
            nbrlist.append(np.array([-1, 0, 1]))
            nbrlist.append(np.array([ 0, 1, 1]))
            nbrlist.append(np.array([ 0,-1, 1]))
            # NNN
            nbrlist.append(np.array([ 2, 0, 0]))
            nbrlist.append(np.array([-2, 0, 0]))
            nbrlist.append(np.array([ 0, 2, 0]))
            nbrlist.append(np.array([ 0,-2, 0]))
            nbrlist.append(np.array([ 0, 0, 2]))
            nbrlist.append(np.array([ 0, 0,-2]))
        else:
            raise ValueError(f'Chosen {self.latt_type} lattice. Currently only FCC lattice is supported.')

        self.nbrlist = nbrlist


    def make_lattice(self, xyz, box):
        """
        Set up site lables on FCC lattice
        """
    
        # create simulation box (fill with -1 denoting sites on on the FCC lattice
        latt = -1*np.ones((tuple(box)), dtype=int)

        # cycle through lattice sites and mark FCC sites with 0 (empty site)
        for r in product(range(box[0]), range(box[1]), range(box[2])):
            if sum(r) % 2 == 0:
                latt[tuple(r)] = 0
    
        # fill lattice sites with atom ids
        for i, r in enumerate(xyz, start=1):
            assert sum(r) % 2 == 0, f'{r} Atom not on FCC lattice!'
            latt[tuple(r)] = i
    
        self.latt = latt
        self.box = box
        self.xyz = xyz
        self.nat = len(self.xyz)

        # Set grain number of each substrate atom to 0
        self.grain = [0 for _ in range(self.nat)]


    def find_neighbors(self, ri):

        # search nearest neighbors
        # to determine stable sites (needs at least 3)
        neighbors = []
        grain_numbers = []
        vacancies = []
        for dr in self.nbrlist[0:12]:
            rj = (ri + dr) % self.box
            iatom = self.latt[tuple(rj)]
            if iatom > 0:
                neighbors.append((iatom, rj))
                grain_numbers.append(self.grain[iatom-1])
            else:
                vacancies.append(rj)

        return neighbors, grain_numbers, vacancies

    def get_grain(self, grain_numbers):

        grain_counts = Counter(grain_numbers)
        if 0 in grain_counts:
            del grain_counts[0]
        
        if any(grain_counts):
            #g_number = sorted(grain_counts, key=grain_counts.__getitem__)[-1]
            g_number = grain_counts.most_common(1)[0][0]
        else:
            g_number = max(self.grain) + 1

        return g_number


    def init_events(self, rates):
 
        print('rates', rates, len(self.xyz), self.box)
        event_list = []

        # Deposition event
        # find sites available for deposition
        # contains at least three neighbors underneath

        # explore or surface sites in the x-y plane
        for ix, iy in product(range(self.box[0]), range(self.box[1])):

            # find z position
            for iz in range(self.box[2]):
                if self.latt[ix, iy, iz] == 0:
                    ri = np.array([ix, iy, iz], dtype=int)
 
                    # count number of neighbors in the target position
                    neighbors, _, _ = self.find_neighbors(ri)

                    # if 3 or more nearest neighbors create a deposition event
                    if len(neighbors) > 2:
                        event = {'type':0, 'rate': rates[0]}
                        event['atom'] = None
                        event['initial'] = None
                        event['final'] = [ix, iy, iz]

                        event_list.append(event)

                    break

        # diffusion events
        for i, ri in enumerate(self.xyz, start=1):
            if ri[2] < 2: continue

            # cycle over neighbor sites
            for dr in self.nbrlist:
                if dr[2] > 0: continue

                rj = (ri + dr) % self.box

                if self.latt[tuple(rj)] == 0:

                    # count number of neighbors in the target position
                    neighbors, _, _ = self.find_neighbors(rj)

                    # if 3 or more nearest neighbors create a diffusion event
                    if len(neighbors) > 2:
                        event = {}
                        event['atom'] = i
                        event['initial'] = ri
                        event['final'] = rj
                        event['type'] = 1
                        event['rate'] = rates[1]
                        event_list.append(event)

        self.event_list =  event_list

        # Initiate event data structures
        self.etree = EventTree(rates, event_list)

        print('events',len(event_list), np.random.choice(event_list, 5))
        print(self.etree.event_tree[-1])
        print('-------')
        print([(e['type'], e['final']) for e in event_list])
        print('-------')


    def move(self, event):
        """
        Perform kmc move given by an event (deposition or diffusion)
        and update event list
        """

        # double check if there are some free spaces (just in case - should
        # follow from zero remaining events)
        if len(self.xyz) == self.box[0]*self.box[1]*self.box[2]/2:
            raise ValueError(f'Lattice is full of atoms, no more events possible.')

        #print('event counts', Counter([e['type'] for e in self.etree.events]))

        old_events = []
        new_events = []

        # deposition event
        if event['type'] == 0:
            ri = event['final']

            # search neighbors and grain numbers
            neighbors, grain_numbers, vacancies = self.find_neighbors(ri)

            # create a new atom
            self.xyz.append(np.array((ix, iy, iz))) 

            # put it on a lattice
            self.latt[tuple(self.xyz[-1])] = len(self.xyz)  # atom number to lattice

            # assign a new grain number to the atom
            g_number = get_grain(grain_numbers)

            # grain number for each atom
            self.grain.append(g_number)

            # Remove the present event

            # Create events around the new atom

            # Adjust events for its neighbors 
            for neighbor in neighbors:
                pass


        else: # diffusion
            ri = event['initial']
            rj = event['final']
            iatom = event['atom']

            # search neighbors of the initial state
            neighbors_old, _, _ = self.find_neighbors(ri)

            # move atom to the new position
            self.latt[tuple(ri)] = 0
            self.latt[tuple(rj)] = iatom+1
            self.xyz[iatom] = rj

            # search neighbors and grain numbers for final state 
            neighbors_new, grain_numbers, vacancies = self.find_neighbors(rj)

            # assign a new grain number to the atom
            g_number = get_grain(grain_numbers)
            
            # grain number for each atom
            self.grain[iatom] = g_number

            # Remove the present event
            old_events.append(event)

            # Change events around the moved atom

            # Add a vacancy-target for old neighbors in the place of ri
            for neighbor in neighbors_old:
                pass

            # Add a vacancy-target for new neighbors in the place of rj
            for neighbor in neighbors_new:
                pass


        # update atom count
        self.nat = len(self.xyz)


        return old_events, new_events

    def get_conf(self):
        return self.xyz, self.box, self.grain

    def advance_time(self):
        """
        Time step of the last event and reset total rates

        Returns
        -------
        dt : float
             Time of the latest event
        """

        dt = -np.log(np.random.random())/self.Rs

        # Save current overall rate (stored in the last node of event tree)
        self.Rs = self.etree.event_tree[-1]

        return dt


    def step(self):
        """
        Perform a KMC step.
        """

        # return a random event
        event = self.etree.find_event()

        # perform a step prescribed by the event and return lists of affected events
        old_events, new_events = self.move(event)
 
        # update binary search tree
        self.etree.update_events(old_events, new_events)
 
