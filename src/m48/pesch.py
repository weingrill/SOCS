#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 8, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
import pylab as plt
import config

from table import Table  # @UnresolvedImport

class Pesch(Table):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        #super(ChildB, self).__init__()
        Table.__init__(self, key = 'No', columns = ['No', 'Ref', 'V', 'B-V'])
        #super(Record, self).__init__()
        #self.key = 'No' 
        self.columns = ['No', 'Ref', 'V', 'B-V']
        #self.table= []
        self._load()
    
    def _load(self):
        filename = config.datapath+'reftable.txt'
        data = np.loadtxt(filename)
        #print data
        for d in data:
            rec = {'No': int(d[0]),
                   'Ref': int(d[1]),
                   'V': float(d[2]),
                   'B-V': float(d[3])}
            self.append(rec)
        
if __name__ == '__main__':
    p = Pesch()
    p.keys()
    print len(p)
    print p
    print p['B-V']
    print set(p['Ref'])
    ref48 = p(Ref=48)
    ref1601 = p(Ref=1601)
    #plt.scatter(p['B-V'], p['V'], edgecolor='none', c=p['Ref'])
    plt.scatter(ref48['B-V'], ref48['V'], edgecolor='none', c='g')
    plt.scatter(ref1601['B-V'], ref1601['V'], edgecolor='none', c='r')
    plt.ylim(plt.ylim()[::-1])
    plt.show()
    