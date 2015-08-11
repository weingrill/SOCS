#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 11, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Stetson(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.data = []
        from datasource import DataSource
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        
    def fromfile(self, filename=None):
        import pyfits
        
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        self.keys = []
        for col in hdulist[1].columns:
            self.keys.append(col.name)
        hdulist.close()
        print self.keys
    
    def todatabase(self):
        for d in self.data:
            print d['ID_1']
            query = """INSERT INTO ngc6633ref (starid, ra, dec, dx, dy, x, y, bmag, bsigma, nb, vmag, vsigma, nv, imag, isigma, ni, coord) 
VALUES ('%(ID_1)s', %(RA)f, %(DEC)f, %(dX)f, %(dY)f, %(X)f, %(Y)f, %(B)f,%(sigmaB)f,%(NB)f,%(V)f,%(sigmaV)f,%(NV)f,%(I)f,%(sigmaI)f,%(NI)f, point(%(RA)f,%(DEC)f));""" % d
            self.wifsip.execute(query, commit=False)
            
        self.wifsip.commit()        

stet = Stetson()
stet.fromfile('/work2/jwe/SOCS/NGC6633/data/NGC6633.fit')
stet.todatabase()        
        