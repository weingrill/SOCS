#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 5, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class StackedPhot(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.data = []
        from datasource import DataSource
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        
    def fromfile(self, filename=None):
        import pyfits
        
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        self.keys = []
        for col in hdulist[1].columns:
            self.keys.append(col.name)
        hdulist.close()
    
    def todatabase(self):
        for d in self.data:
            print d['starid']
            record = {}
            for key in self.keys:
                record[key] = d[key]
            query = """UPDATE ngc2236 SET vmag = %(VMAG)f, bmag = %(BMAG)f, bv = %(BMAG)f-%(VMAG)f, nv = %(V_FLAGS)d, nb = %(B_FLAGS)d, vmag_err = %(VMAG_ERR)f, bmag_err = %(BMAG_ERR)f, member=TRUE WHERE starid = '%(starid)s';""" % record
            self.wifsip.execute(query, commit=False)
            
        self.wifsip.commit()        

vphot = StackedPhot()
vphot.fromfile('/work2/jwe/NGC2236/ngc2236stackedphot.fits')
vphot.todatabase()        