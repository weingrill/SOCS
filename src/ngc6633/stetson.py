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
    
    def calibrate(self):
        import numpy as np
        ffield = {'V': 'vmag', 'B': 'bmag', 'I': 'imag'}
        for filtercol in ['V','B','I']:
            query = """SELECT objid 
            FROM frames 
            WHERE object LIKE 'NGC 6633 BVI %%' 
            AND filter = '%s'
            AND corr >=0;""" % filtercol
            frames = self.wifsip.query(query)
            for frame in frames:
                print frame[0]
                params = {'objid': frame[0],
                          'filterfield': ffield[filtercol]}
                query = """SELECT starid, mag_isocor, magerr_isocor, %(filterfield)s, mag_isocor - %(filterfield)s
                FROM phot, ngc6633ref
                WHERE objid = '%(objid)s'
                AND circle(phot.coord,0) <@ circle(ngc6633ref.coord, 0.3/3600.0)
                AND NOT %(filterfield)s IS NULL
                AND flags<4;""" % params
                result = self.wifsip.query(query)
                if len(result)>0:
                    ocs = np.array([r[4] for r in result])
                    errs = np.array([r[2] for r in result])
                    weights=1/errs
                    for r in result:
                        print r
                    corr = np.average(ocs, weights=weights)
                    print 'corr = %.4f' % corr
                    params['corr'] = corr
                    query = """UPDATE frames SET corr=%(corr)f WHERE objid ='%(objid)s';""" % params
                    self.wifsip.execute(query)

stet = Stetson()
#stet.fromfile('/work2/jwe/SOCS/NGC6633/data/NGC6633.fit')
#stet.todatabase()        
stet.calibrate()