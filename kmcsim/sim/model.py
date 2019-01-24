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
from collections import Counter, defaultdict
#import events
from .events import EventTree

class KMCModel:
    """Class managing kmc moves and event modifications"""

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
        Set up site labels on FCC lattice
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
        for dr in self.nbrlist[0:12]:
            rj = tuple((ri + dr) % self.box)
            iatom = self.latt[rj]
            neighbors.append((iatom, rj))
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

    def find_events(self, j, rj):
        """
        Finds events for site j at rj.
        Should be used in init_events
        """

        ix, iy, iz = rj

        events_found = []
        iatom = self.latt[ix, iy, iz]

        # vacancy, test for possibility of deposition event
        if iatom == 0:
           # rj = np.array(rj, dtype=int)

            # count number of atomic neighbors in the target position
            _, grain_numbers_j = self.find_neighbors(rj)

            # if 3 or more nearest neighbors with grain IDs, create a deposition event
            if len(grain_numbers_j) > 2:
                events_found.append((0, rj, rj, rj))

        # atom, find diffusion events
        elif iatom > 0:

            # search for possible diffusion events

            # explore neighborhood of atom j
            neighbors_j, grain_numbers_j = self.find_neighbors(rj)

            for k, rk in neighbors_j:

                # find if vacancy is a good destination spot
                if k == 0:
                    _, grain_numbers_k = self.find_neighbors(rk)

                    # if number of real atoms 3 or more, make vacancy available as
                    # a destination for deposition and diffusion
                    if len(grain_numbers_k)-1 > 2:
                        # do not diffuse upward
                        if rk[2] > rj[2]:
                            continue
                        events_found.append((1, j, rj, rk))

        return events_found


    def init_events(self, rates):
 
        rates = np.array(rates)

        # structure to store event information
        event_list = [[] for _ in range(rates.shape[0])]

        # dictionary to store references to event_list
        site_dict = defaultdict(list)

        # Deposition event - find vacant sites available for deposition
        for ix, iy in product(range(self.box[0]), range(self.box[1])):

            # find z position
            for iz in range(self.box[2]):
                if self.latt[ix, iy, iz] == 0:
                    ri = np.array([ix, iy, iz], dtype=int)
 
                    # count number of neighbors in the target position
                    neighbors, grain_numbers = self.find_neighbors(ri)

                    # if 3 or more nearest neighbors with grain IDs, create a deposition event
                    if len(grain_numbers) > 2:
                        event = {'type':0}
                        event['atom'] = 0
                        event['initial'] = (ix, iy, iz)
                        event['final'] = (ix, iy, iz)
                        event_list[0].append(event)

                        # add event information to the site
                        site_dict[(ix, iy, iz)].append((0, len(event_list[0])-1))
                    break

        # diffusion events for actual atoms (i.e., atom id > 0)
        for i, ri in enumerate(self.xyz, start=1):
            if ri[2] < 2: continue

            # cycle over neighbor sites
            for dr in self.nbrlist:

                # do not diffuse upward
                if dr[2] > 0: continue

                rj = (ri + dr) % self.box

                # is vacancy in the neighborhood of atom i?
                if self.latt[tuple(rj)] == 0:

                    # explore neighborhood of the target vacancy
                    neighbors, grain_numbers = self.find_neighbors(rj)

                    # if 3 or more nearest neighbors with grain IDs present, create a diffusion event
                    if len(grain_numbers) > 2:
                        event = {'type':1}
                        event['atom'] = 1
                        event['initial'] = ri
                        event['final'] = rj
                        event_list[1].append(event)

                        # add event information to the site
                        site_dict[i].append((1, len(event_list[1])-1))

        self.event_list =  event_list
        self.site_dict = site_dict

        # Get dictionary of event type counts
        n_events = np.array([len(e) for e in self.event_list])
        print('Number of events:', n_events)

        # Initiate event data structures
        self.etree = EventTree(rates)
        self.etree.update_events(n_events)

        #print('events',len(event_list), np.random.choice(np.random.choice(event_list, 2), 5))
        #print(self.etree.event_tree[-1])
        #print('-------')
        #print([[(e['type'], e['final']) for e in et] for et in event_list])
        #print('-------')
        #print('site dict', site_dict)


    def move(self, event_type, event_number):
        """
        Perform kmc move given by an event (deposition or diffusion)
        and update event list
        """

        # double check if there are some free spaces (just in case - should
        # follow from zero remaining events)
        if len(self.xyz) == self.box[0]*self.box[1]*self.box[2]/2:
            raise ValueError(f'Lattice is full of atoms, no more events possible.')

        event = self.event_list[event_type][event_number]
        old_events = []
        new_events = []
        n_events = self.etree.n_events

        print('# event types:', event_type, [len(el) for el in self.event_list])
        print('atom numbers:', len(self.xyz), len(set(self.grain)))
        print('last atom coords:', self.xyz[-1])

        for i in range(len(self.event_list)):
            #print('# start events: #', i, len(self.event_list[i]), n_events[i])
            assert len(self.event_list[i]) == n_events[i], f'Number of events of type {i} does not match: {len(self.event_list[i])} vs. {n_events[i]}' 

        # deposition event
        if event_type == 0:
            ri = event['final']

            # create a new atom
            self.xyz.append(np.array(ri)) 
            iatom = len(self.xyz)

            # put it on a lattice
            self.latt[tuple(self.xyz[-1])] = iatom # atom number to lattice

            # search neighbors and grain numbers
            neighbors, grain_numbers = self.find_neighbors(ri)

            # assign a new grain number to the atom
            self.grain.append(self.get_grain(grain_numbers))

            # Identify old events for removal

            # remove the current deposition event
            old_events.append((event_type, event_number))
            # ... and the associated dictionary of site events
            #!!! FIX
            #del self.site_dict[ri]
            n_events[0] -= 1

            # find events of the moved atom
            events_found = self.find_events(iatom, ri)

            # update site dict and new_events list cycle through new events
            for ev in events_found:
                new_events.append(ev)
                n_events[ev[0]] += 1


            # remove all old events of the old and new neighbors
            # and add new events
            for j, rj in neighbors:

                # remove all current events of neighbor j
                for et, en in self.site_dict[j]:
                    old_events.append((et, en))
                    self.site_dict[j].remove((et, en))
                    n_events[et] -= 1

                if len(self.site_dict[j]) == 0:
                    del self.site_dict[j]

                events_found = self.find_events(j, rj)

                # update site dict and new_events list cycle through new events
                for ev in events_found:
                    new_events.append(ev)
                    n_events[ev[0]] += 1

        elif event_type == 1: # diffusion
            r0 = event['initial']
            ri = event['final']

            # identify atom (to access associated events)
            iatom = self.latt[r0]

            # search neighbors of the initial state
            neighbors_old, _ = self.find_neighbors(r0)

            # create move atom to the new position
            self.latt[tuple(r0)] = 0
            self.latt[tuple(ri)] = iatom
            print('diff', r0, ri, iatom, self.xyz[iatom-1])
            self.xyz[iatom-1] = ri

            # remove all current events of atom iatom
            for et, en in self.site_dict[iatom]:
                self.site_dict[iatom].remove((et, en))
                old_events.append((et, en))
                n_events[et] -= 1

            if len(self.site_dict[iatom]) == 0:
                del self.site_dict[iatom]

            # remove all current events of the destination vacancy
            for et, en in self.site_dict[ri]:
                self.site_dict[ri].remove((et, en))
                old_events.append((et, en))
                n_events[et] -= 1

            # find events of the moved atom
            events_found = self.find_events(iatom, ri)

            # update site dict and new_events list cycle through new events
            for ev in events_found:
                new_events.append(ev)
                n_events[ev[0]] += 1

            # search neighbors and grain numbers for final state 
            neighbors_new, grain_numbers = self.find_neighbors(ri)

            # assign a new grain number to the atom
            self.grain[iatom-1] = self.get_grain(grain_numbers)

            # combine all neighbors (old and new)
            neighbors = list(set(neighbors_old + neighbors_new))

            # remove all old events of the old and new neighbors
            # and add new events
            for j, rj in neighbors:

                # remove all current events of neighbor j
                for et, en in self.site_dict[j]:
                    self.site_dict[j].remove((et, en))
                    old_events.append((et, en))
                    n_events[et] -= 1

                if len(self.site_dict[j]) == 0:
                    del self.site_dict[j]

                events_found = self.find_events(j, rj)

                # update site dict and new_events list cycle through new events
                for ev in events_found:
                    new_events.append(ev)
                    n_events[ev[0]] += 1

        # update events lists with old_events and new_events
        print('# end events:', event_type, len(old_events), len(new_events))

        # identify indices to be deleted
        delete_events = [[] for _ in range(len(self.event_list))]
        for ev in old_events:
            delete_events[ev[0]].append(ev[1])

        for i in range(len(self.event_list)):
            self.event_list[i] = [ev for k, ev in enumerate(self.event_list[i]) if k not in delete_events[i]]

        for ev in new_events:
            n_event = {'type':ev[0]}
            n_event['atom'] = ev[1]
            n_event['initial'] = ev[2]
            n_event['final'] = ev[3]
            self.event_list[ev[0]].append(n_event)
            self.site_dict[ev[1]].append((ev[0], len(self.event_list[ev[0]])-1 ))

        for i in range(len(self.event_list)):
            assert len(self.event_list[i]) == n_events[i], f'Number of events of type {i} does not match: {len(self.event_list[i])} vs. {n_events[i]}' 

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
        print('ee', event_type, event_number)

        # perform a step prescribed by the event and return lists of affected events
        n_events = self.move(event_type, event_number)

        # update binary search tree
        self.etree.update_events(n_events)
 
