#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 25, 2016

@author: Joerg Weingrill <jweingrill@aip.de>

Determine the Halpha index of M67 stars
'''
import numpy as np
import config

class M67Alpha(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self._loaddata()
        
    def _loaddata(self, filename = config.datapath+'M67_by_alpha.tst'):
        '''
        loads the datafile with matched ymag, bmag, hanmag, hawmag and errors
        '''
        
        colnames = ['NUMBER_1','bmag','bmag_err','ra','dec','ymag','ymag_err','hawmag','hawmag_err','hanmag','hanmag_err','Index']
        
        self.data = np.genfromtxt(filename, dtype=None, delimiter='\t', skip_header=24, skip_footer=1, names=colnames)
        print self.data['NUMBER_1']
        self.bymag = self.data['bmag'] - self.data['ymag'] + 0.63830
        self.alpha = (self.data['hanmag'] - self.data['hawmag'])*2 + 1.9
        
        
    def plot(self):
        '''
        plot the data
        '''
        import matplotlib.pyplot as plt
        plt.scatter(self.bymag, self.alpha)
        plt.axhline(2.06, linestyle='--', label='zero EW')
        plt.xlabel('b - y')
        plt.ylabel('$\\alpha$ [mag]')
        plt.xlim(0.2, 1.2)
        plt.ylim(2.5,1.3)
        plt.show()
    
    
if __name__ == '__main__':
    m = M67Alpha()
    m.plot()