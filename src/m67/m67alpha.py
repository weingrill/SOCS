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
        self._loadmembers()
        self._loadgood()
        
    def _loadtstfile(self, filename):
        '''
        loads a TST File written with STIL v3.0-9
        '''
        with open(filename, 'rt') as f:
            lines = f.readlines()
        #print lines[22]
        for i,line in enumerate(lines):
            if line[0:10] == '# This TST':
                break
        colnames = lines[i+2].rstrip('\n').split('\t')
        return np.genfromtxt(filename, dtype=None, delimiter='\t', skip_header=i+4, skip_footer=1, names=colnames)
        
        
    def _loaddata(self, filename = config.datapath+'M67_by_alpha.tst'):
        '''
        loads the datafile with matched ymag, bmag, hanmag, hawmag and errors
        '''
        self.data = self._loadtstfile(filename)
        self.bymag = self.data['bmag'] - self.data['ymag'] + 0.63830
        self.alpha = (self.data['hanmag'] - self.data['hawmag'])*2 + 1.9

    def _loadmembers(self, filename = config.datapath+'M67_by_alpha_K2_members.tst'):
        '''
        loads the datafile with matched ymag, bmag, hanmag, hawmag and errors
        '''
        self.members = self._loadtstfile(filename)
        #self.bymag = self.data['bmag'] - self.data['ymag'] + 0.63830
        #self.alpha = (self.data['hanmag'] - self.data['hawmag'])*2 + 1.9

    def _loadgood(self, filename = config.datapath+'M67_by_alpha_K2_good.tst'):
        '''
        loads the datafile with matched ymag, bmag, hanmag, hawmag and errors
        '''
        self.good = self._loadtstfile(filename)
        #self.bymag = self.data['bmag'] - self.data['ymag'] + 0.63830
        #self.alpha = (self.data['hanmag'] - self.data['hawmag'])*2 + 1.9
        
        
    def plot(self, filename=None):
        '''
        plot the data
        '''
        import matplotlib.pyplot as plt
        
        goodby = self.good['bmag']-self.good['ymag'] + 0.63830
        goodalpha = (self.good['hanmag'] - self.good['hawmag'])*2 + 1.9
        
        memby = self.members['bmag']-self.members['ymag'] + 0.63830
        memalpha = (self.members['hanmag'] - self.members['hawmag'])*2 + 1.9
        
        plt.figure(figsize=[6,6])
        plt.scatter(self.bymag, self.alpha, color='k', alpha=0.5, edgecolor='none',label='field')
        plt.scatter(memby, memalpha, color='b', label='members')
        plt.scatter(goodby, goodalpha, color='r', label='rotators')
        
        plt.plot( 0.4105, 2.196, ',', marker='$\\bigodot$', markersize=8, label='sun', color='orange')
        
        plt.axhline(2.06, linestyle='--', label='zero EW', color='g')
        plt.xlabel('b - y')
        plt.ylabel('$\\alpha$ [mag]')
        #plt.xlim(0.2, 1.2)
        plt.xlim(0.3, 0.8)
        plt.ylim(2.5,1.3)
        plt.minorticks_on()
        plt.legend(scatterpoints=1, fontsize='smaller', numpoints=1)
        plt.tight_layout()
        if type(filename) is str:
            plt.savefig(filename)
        elif filename is None:
            plt.show()
        else: 
            raise(TypeError)
        plt.close()
    
    
if __name__ == '__main__':
    m = M67Alpha()
    m.plot(config.resultpath+'M67_alpha.pdf')