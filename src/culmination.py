#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2018 AIP"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__="2018-06-18"

import numpy as np
import matplotlib.pyplot as plt 
 
class Culmination(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.ra = np.array([84.59, 282.84, 29.28, 132.85, 105.69])
        self.month = np.array([12.5, 7.1, 10.83, 1.9, 1.13])
    
    def plot(self):
        month = ((np.arange(24) - 5.86666 +24.) % 24) / 2.
        plt.plot(self.ra / 15., self.month, 'ko')
        plt.plot(month)
        plt.xlim(24.0,0.0)
        plt.xticks(np.arange(24))
        plt.yticks(np.arange(13))
        plt.show()
            
            
c = Culmination()
c.plot()