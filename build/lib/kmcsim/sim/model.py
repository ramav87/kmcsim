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
import random
import numpy as np
from itertools import product
from collections import Counter, defaultdict
#import events
from .events import EventTree

class KMCModel:
    """Class managing kmc moves and event modifications"""

    def __init__(self, latt_type):
        self.latt_type = latt_type
        self.__setup_neighbors()
        self.verbose=False

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
            raise ValueError('Chosen {} lattice. Currently only FCC lattice is supported.'.format(self.latt_type))

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
            assert sum(r) % 2 == 0, '{} Atom not on FCC lattice!'.format(r)
            latt[tuple(r)] = i
    
        self.latt = latt
        self.box = box
        self.xyz = xyz
        self.nat = len(self.xyz)

        # Set grain number of each substrate atom to 0
        self.grain = [0 for _ in range(self.nat)]


    def find_neighbors(self, ri):
        """
        Find neighbors and return ids usable in site_dict
        """

        # search nearest neighbors
        # to determine stable sites (needs at least 3)
        neighbors = []
        grain_numbers = []
        for dr in self.nbrlist[0:12]:
            rj = tuple((np.array(ri) + dr) % self.box)
            iatom = self.latt[rj]
            neighbors.append(rj)
            if iatom > 0:
                grain_numbers.append(self.grain[iatom-1])

        return neighbors, grain_numbers

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

    def find_events(self, rj):
        """
        Finds events for site j at rj.
        Should be used in init_events
        """

        ix, iy, iz = rj
        iatom = self.latt[ix, iy, iz]
        #print('rj', type(rj), rj, iatom)

        events_found = []

        # vacancy, test for possibility of deposition event
        if iatom == 0:
           # rj = np.array(rj, dtype=int)

            # count number of atomic neighbors in the target position
            _, grain_numbers_j = self.find_neighbors(rj)

            # if 3 or more nearest neighbors with grain IDs, create a deposition event
            if len(grain_numbers_j) > 2:
                events_found.append((0, ix, iy, iz, ix, iy, iz))

        # atom, find diffusion events
        elif iatom > 0:

            # search for possible diffusion events

            # explore neighborhood of atom j
            neighbors_j, grain_numbers_j = self.find_neighbors(rj)

            for rk in neighbors_j:

                # find if vacancy is a good destination spot
                if self.latt[rk] == 0:
                    _, grain_numbers_k = self.find_neighbors(rk)

                    # if number of real atoms 3 or more, make vacancy available as
                    # a destination for deposition and diffusion
                    if len(grain_numbers_k)-1 > 2:

                        if len(np.unique(grain_numbers_k)) <=2:
                            # do not diffuse upward
                            if rk[2] > rj[2]:
                                continue
                            events_found.append((1, ix, iy, iz, rk[0], rk[1], rk[2]))
                        else:
                            # do not diffuse upward
                            if rk[2] > rj[2]:
                                continue
                            events_found.append((2, ix, iy, iz, rk[0], rk[1], rk[2]))


        return events_found


    def init_events(self, rates):
 
        rates = np.array(rates)

        # structure to store event information list of sets, containing
        # events dictionary keys
        event_list = [set() for _ in range(rates.shape[0])]

        # dictionary to store references to event_list
        site_dict = defaultdict(list)

        # Deposition event - find vacant sites available for deposition
        for ix, iy in product(range(self.box[0]), range(self.box[1])):

            # find z position
            for iz in range(self.box[2]):
                if self.latt[ix, iy, iz] == 0:
                    t_ri = (ix, iy, iz)
 
                    # count number of neighbors in the target position
                    neighbors, grain_numbers = self.find_neighbors(t_ri)

                    # if 3 or more nearest neighbors with grain IDs, create a deposition event
                    if len(grain_numbers) > 2:
                        event_tuple = (0, ix, iy, iz, ix, iy, iz)
                        event_list[0].add(event_tuple)

                        # add event information to the site
                        site_dict[t_ri].append(event_tuple)
                    break

        # diffusion events for actual atoms (i.e., atom id > 0)
        for i, ri in enumerate(self.xyz, start=1):
            if ri[2] < 2: continue

            # cycle over neighbor sites
            for dr in self.nbrlist:

                # do not diffuse upward
                if dr[2] > 0: continue

                rj = tuple((ri + dr) % self.box)

                # is vacancy in the neighborhood of atom i?
                if self.latt[rj] == 0:

                    # explore neighborhood of the target vacancy
                    neighbors, grain_numbers = self.find_neighbors(rj)

                    # if 3 or more nearest neighbors with grain IDs present, create a diffusion event
                    if len(grain_numbers) > 2:
                        # Add cross-grain diffusion. This is true if the grains themselves are different.
                        if len(np.unique(grain_numbers)) > 2:
                            event_tuple = (2, ri[0], ri[1], ri[2], rj[0], rj[1], rj[2])
                            event_list[2].add(event_tuple)
                            # add event information to the site
                            site_dict[tuple(ri)].append(event_tuple)
                        else:
                            event_tuple = (1, ri[0], ri[1], ri[2], rj[0], rj[1], rj[2])
                            event_list[1].add(event_tuple)
                            # add event information to the site
                            site_dict[tuple(ri)].append(event_tuple)



        self.event_list =  event_list
        self.site_dict = site_dict

        # Get dictionary of event type counts
        n_events = np.array([len(e) for e in self.event_list])
        print('Number of events:', n_events)

        # Initiate event data structures
        self.etree = EventTree(rates)
        self.etree.update_events(n_events)


    def move(self, event_type, event_number):
        """
        Perform kmc move given by an event (deposition or diffusion)
        and update event list
        """

        # double check if there are some free spaces (just in case - should
        # follow from zero remaining events)
        if len(self.xyz) == self.box[0]*self.box[1]*self.box[2]/2:
            raise ValueError('Lattice is full of atoms, no more events possible.')

        # find a tuple containing information about the selected event
        event = tuple(self.event_list[event_type])[event_number]
        old_events = []
        new_events = []
        n_events = self.etree.n_events

        #print('# event:', event, 'ev#', [len(el) for el in self.event_list], end='')
        #print('at#',len(self.xyz), 'gr#', len(set(self.grain)), 'lxyz', self.xyz[-1])

        for i in range(len(self.event_list)):
            assert len(self.event_list[i]) == n_events[i], 'Start: Number of events of type {} does not match: {} vs. {}'.format(i,len(self.event_list[i]),n_events[i])

        # deposition event
        if event_type == 0:
            t_ri = event[4:7]
            ri = np.array(t_ri)

            # create a new atom
            self.xyz.append(ri) 
            iatom = len(self.xyz)

            # put it on a lattice
            self.latt[t_ri] = iatom # id for the site properties with atom id and list of events

            # search neighbors and grain numbers
            neighbors, grain_numbers = self.find_neighbors(t_ri)

            # assign a new grain number to the atom
            self.grain.append(self.get_grain(grain_numbers))

            # Identify old events for removal
            # remove the current deposition event
            old_events.extend(self.site_dict[t_ri])
            # ... and the associated dictionary of site events
            del self.site_dict[t_ri]

            # find diffusion events of the deposited atom
            events_found = self.find_events(t_ri)
            new_events.extend(events_found)

            # remove all old events of the new neighbors and add their new events
            for t_rj in neighbors:

                if t_rj == t_ri:
                    continue

                # remove all current events of neighbor j
                old_events.extend(self.site_dict[t_rj])
                del self.site_dict[t_rj]

                # add new events of neighbor j
                events_found = self.find_events(t_rj)
                new_events.extend(events_found)


        else: # diffusion
            t_r0 = event[1:4] # initial position
            t_ri = event[4:7] # final position
            r0 = np.array(t_r0)
            ri = np.array(t_ri)

            # identify atom (to access associated events)
            iatom = self.latt[t_r0]

            # remove all current events of atom iatom
            old_events.extend(self.site_dict[t_r0])
            del self.site_dict[t_r0]

            # remove all current events of the destination vacancy
            old_events.extend(self.site_dict[t_ri])
            del self.site_dict[t_ri]

            # search neighbors of the initial state
            neighbors_old, _ = self.find_neighbors(t_r0)

            # move atom to the new position
            self.latt[t_r0] = 0
            self.latt[t_ri] = iatom
            self.xyz[iatom-1] = ri

            # find events of the moved atom
            events_found = self.find_events(ri)

            # update site dict and new_events list cycle through new events
            new_events.extend(events_found)

            # search neighbors and grain numbers for final state 
            neighbors_new, grain_numbers = self.find_neighbors(t_ri)

            # assign a new grain number to the atom
            self.grain[iatom-1] = self.get_grain(grain_numbers)

            # remove all old events of the old and new neighbors
            # and add new events
            for t_rj in set(neighbors_old + neighbors_new):

                if t_rj == t_ri or t_rj == t_r0:
                    continue

                old_events.extend(self.site_dict[t_rj])
                del self.site_dict[t_rj]

            # add new events of neighbor j
            events_found = self.find_events(t_rj)
            new_events.extend(events_found)


        # update events lists with old_events and new_events
        o_events = 0
        for i in range(len(self.event_list)):
            o_events += len(self.event_list[i])

        for ev in old_events:
            self.event_list[ev[0]].remove(ev)

        for ev in new_events:
            if ev in self.event_list[ev[0]]:
                print('present', ev)
            self.event_list[ev[0]].add(ev)
            self.site_dict[ev[1:4]].append(ev)

        n_events = []
        for i in range(len(self.event_list)):
            n_events.append(len(self.event_list[i]))

        df = len(new_events) - len(old_events)
        sm = sum(n_events) - o_events

        print('Num new events: {}, num Old events: {}'.format(len(new_events), len(old_events)))
        print('Num of events: {}'.format(sum(n_events)))
        
        assert sm == df, "Number of new-old events does not match: {0} {1}".format(sm, df)

        # update atom count
        self.nat = len(self.xyz)

        return n_events


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

        dt = -np.log(np.random.random())/self.etree.Rs

        return dt


    def step(self):
        """
        Perform a KMC step.
        """

        # return a random event (based on their frequency)
        event_type, event_number = self.etree.find_event()

        # perform a step prescribed by the event and return lists of affected events
        n_events = self.move(event_type, event_number)

        # update binary search tree
        self.etree.update_events(n_events)
 
