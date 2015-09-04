#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 5, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config

class Yadav(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.data = []
        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        
    def fromfile(self, filename=None):
        import pyfits
        
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        self.keys = []
        for col in hdulist[1].columns:
            self.keys.append(col.name)
        hdulist.close()
    
    def todatabase(self):
        import numpy as np
        for d in self.data:
            print d['seq']
            record = {}
            for key in self.keys:
                record[key] = d[key]
                if np.isnan(record[key]):
                    record[key] = 'NULL'
            query = """INSERT INTO m67 (seq, ra, dec, vmag, bmag, icmag, pmra, pmdec, pmb, hrv) 
VALUES (%(Seq)d, %(RAJ2000)f, %(DEJ2000)f, %(Vmag)s, %(Bmag)s, %(Icmag)s, %(pmRA)f, %(pmDE)f, %(Pmb)f, %(HRV)s);""" % record
            self.wifsip.execute(query, commit=False)
            
        self.wifsip.commit()        

y = Yadav()
y.fromfile(config.datapath+'/Yadav.fit')
y.todatabase()        
