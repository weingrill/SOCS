#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 8, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
import pylab as plt


class Record(list):
    '''
    simulates data records
    '''
    def __init__(self):
        self.table=[]
    
    def __getitem__(self, key):
        if type(key) is int:
            return self.table[key]
        
        if type(key) is str:
            return [t[key] for t in self.table]

class Pesch(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
        self.table=[]
        self._loadpesch()
    
    def _loadpesch(self):
        filename = '/work2/jwe/M48/data/reftable.txt'
        data = np.loadtxt(filename)
        print data
        for d in data:
            rec = {'No': int(d[0]),
                   'Ref': int(d[1]),
                   'V': float(d[2]),
                   'B-V': float(d[3])}
            self.table.append(rec)
        plt.scatter([t['B-V'] for t in self.table],[t['B-V'] for t in self.table])
        plt.show()
        
if __name__ == '__main__':
    p = Pesch()