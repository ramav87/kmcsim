#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import numpy as np
from collections import Counter

class EventTree:
    """
    Class maintaining a binary tree for random event type lookup
    and arrays for choosing specific event. 
    """

    def __init__(self, rates, events):

        self.rates = rates
        self.events = events
        self.__setup()


    def __build_tree(self, e_ratio):

        self.event_tree = []

        # create event ratio array level 0 - bottom
        if len(e_ratio) % 2 == 1:
            e_ratio.extend([0.0])

        # create the bottom level (rates*numbers)
        self.event_tree.append(np.array(e_ratio))

        # create partial summs (iteratively) up to the 2nd highest level
        while len(e_ratio) > 2:
            e_ratio = [e_ratio[i]+e_ratio[i+1] for i in range(0, len(e_ratio), 2)]
            if len(e_ratio) % 2 == 1:
                e_ratio.extend([0.0])

            self.event_tree.append(np.array(e_ratio))

        # create top level = sum of all rates
        self.event_tree.append(np.array(sum(e_ratio)))
        

    def __setup(self):

        # Get dictionary of event type counts
        e_counts = Counter([e['type'] for e in self.events])
        print(e_counts)

        # create a list of events based on event types
        self.event_counts = [[] for _ in range(len(self.rates))]
        for e in self.events:
            self.event_counts[e['type']].append(e)

        e_ratio = [e_counts.get(t, 0)*r for t, r in enumerate(self.rates)]

        print('e_ratio', e_ratio)
        self.__build_tree(e_ratio)


    def update_events(self, old_events, new_events):
        """
        Update tree: remove old events and add new events
        """
        pass


    def find_event(self):
        """Find and return an event"""

        # generate a random number [0,Rs)
        q = self.Rs*np.random.random()

        # cycle through levels (top->down)
        # start with top-level child (k-2) end with level above bottom (1)
        j = 0
        for k in range(len(self.event_tree)-2, 0, -1):
            # left child value
            left = self.event_tree[k][j]

            if q < left:
                j = 2*j
            else:
                q -= left
                j = 2*j + 1
        
        # bottom level - return selected event type
        if q < self.event_tree[0][j]:
            event_type = self.events[j]
        else:
            event_type = self.events[j+1]

        # select a random event index of a given type 
        event_number = np.random.randint(len(self.event_counts[event_type]))

        # get the event object        
        event = event_counts[event_type][event_number]

        return event

