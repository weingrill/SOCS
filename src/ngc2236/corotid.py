#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 31, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class CoRoTId(object):
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
            filename = '/work2/jwe/NGC2236/NGC2236corotid.fits'
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        hdulist.close()
    
    def todatabase(self):
        for d in self.data:
            print d['starid'],d['corotid_2']
            self.wifsip.execute("""UPDATE ngc2236
                          SET corotid = %d
                          WHERE starid = '%s';""" % (d['corotid_2'], d['starid']), commit=False)
        self.wifsip.commit()        

cid = CoRoTId()
cid.fromfile()
cid.todatabase()