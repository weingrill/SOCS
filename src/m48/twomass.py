#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 31, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class TwoMass(object):
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
        
        if filename is None:
            filename = '/work2/jwe/m48/data/2MASS.fits'
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        hdulist.close()
    
    def todatabase(self):
        for d in self.data:
            print d['starid'],d['ra_cone'],d['dec_cone']
            self.wifsip.execute("""UPDATE m48stars
                          SET ra = %f, dec = %f
                          WHERE starid = '%s';""" % (d['ra_cone'],d['dec_cone'], d['starid']), commit=False)
        self.wifsip.commit()        

tm = TwoMass()
tm.fromfile()
tm.todatabase()