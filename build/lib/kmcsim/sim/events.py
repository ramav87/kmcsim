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

    def __init__(self, rates):

        # array of reaction rates 
        self.rates = np.array(rates)

        # array for number of events of a given type (same length as rates)
        self.n_events = np.zeros(self.rates.shape, dtype=int)

        # setup event selection data structures
        self.__setup_tree()


    def __setup_tree(self):
        """
        Builds an empty binary search tree structure (list of arrays) for random event type selection.
        The values are filled in the update_events method
        """

        # create an event decision tree structure
        self.event_tree = []

        nrates = len(self.rates)

        # create arrays of exponentially decreasing length divisible by 2
        while nrates > 1:

            if nrates % 2 == 1:
                nrates += 1

            self.event_tree.append(np.zeros((nrates), dtype=float))

            nrates //= 2

        # top node - total rate
        self.event_tree.append(np.zeros((1), dtype=float))

        self.kmax = len(self.event_tree)

        # check tree structure
        for i, t in enumerate(self.event_tree):
            pass
            #print('tree level', i, t.shape)


    def update_events(self, n_events):
        """
        Update tree with new values, if needed
        """

        assert len(n_events) == len(self.rates), 'Rates and n_event lists do not match'

        # if no changes, return
        #if np.array_equal(self.n_events, n_events):
        #    return

        self.n_events[:] = n_events

        # fill the base level (leaves)
        self.event_tree[0][:self.rates.shape[0]] = self.rates*self.n_events

        # create partial summs up to the top level
        for k in range(1, self.kmax):
            self.event_tree[k][:] = [self.event_tree[k-1][i] + self.event_tree[k-1][i+1] for i in range(0, self.event_tree[k-1].shape[0], 2)]

        self.Rs = self.event_tree[-1][0]

        #for i, t in enumerate(self.event_tree):
        #    print('tree level', i, t)
       
    def update_global_rate(self, new_rates):
        """Function that updates the global rates, e.g. if we want to change sim parameters"""

        self.rates = new_rates
        #Update the tree with the new global rates
        #Remember to update events if you update rates

    def find_event(self):
        """Find and return an event"""

        #for i, t in enumerate(self.event_tree):
        #    print('tree level runtime', i, t)

        # generate a random number [0,Rs)
        q = self.Rs*np.random.random()

        # cycle through levels (top->down)
        # start with top-level child (k-2) end with level above bottom (1)
        j = 0
        for k in range(self.kmax-2, -1, -1):

            # left child value
            left = self.event_tree[k][j]

            if q < left:
                j = 2*j
            else:
                q -= left
                j = 2*j + 1
        
        event_type = j


        # select a random event index of a given type
        print('Event type: {}, N_Events: {}'.format(event_type, self.n_events))
        event_number = np.random.randint(low= 0, high = self.n_events[event_type])
        #print('event tree:', event_type, event_number)


        return event_type, event_number

