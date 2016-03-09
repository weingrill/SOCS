#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Mar 7, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

from frame import Frame
import numpy as np
from matplotlib import pyplot as plt

def refstars():
    from datasource import DataSource
    
    table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
    query = "SELECT * FROM referencestars order by starid;"
    data = table.query(query)
    
    columns = data[0].keys()
    print data[0].values()
    dtypes = [type(d) for d in data[0].values()]
    print dtypes
    arraydata = []
    for star in data:
        arraydata.append(tuple(star))
    return np.array(arraydata, dtype = zip(columns, dtypes))
    

class ReferenceStars(object):
    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
        columns = ['starid','ra','dec','"b-y"','m1','c1','beta']
        query = "SELECT %s FROM referencestars order by starid;" % ', '.join(columns)
        data = table.query(query)
        columns[3] = 'b-y'
        
        d0 = data[0]
        #for c in columns: print c
        dtypes = ['|S11','f4','f4','f4','f4','f4','f4']
        for c in columns:
            print c
            dtypes.append(type(d0[c]))
        #dtypes = [type(d) for d in data[0].values()]
        
        for c, d in zip(columns, dtypes): print c,d
        arraydata = []
        for star in data:
            arraydata.append(tuple(star))
        self.stars = np.array(arraydata, dtype = zip(columns, dtypes))
    

class CalStromgren(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.loadframe()
        self.loadref()
    
    def loadframe(self):
        self.bframe = Frame('20151126A-0008-0006')    
        self.yframe = Frame('20151126A-0008-0005')    
        
    def loadref(self):
        self.reference = ReferenceStars()
        pass
    
    def match(self):
        from astropy.coordinates import SkyCoord,search_around_sky   # @UnresolvedImport
        from astropy import units as u
        
        # fill in the coordinates and b-y
        nstars = len(self.reference.stars['ra'])
        rra2 = self.reference.stars['ra']
        rdec2 = self.reference.stars['dec']
        refcoord = SkyCoord(ra=rra2*u.degree, dec=rdec2*u.degree)
        cby = self.reference.stars['b-y']
        b = np.ones(nstars) * np.nan
        y = np.ones(nstars) * np.nan
        
        bra1 = self.bframe.stars['alphawin_j2000']
        bdec1 = self.bframe.stars['deltawin_j2000']
        bcoords = SkyCoord(ra=bra1*u.degree, dec=bdec1*u.degree)  
        bfluxes = self.bframe.stars['flux_auto']/self.bframe.expt
        
        bidx, refidx, _, _ = search_around_sky(bcoords, refcoord, 0.6*u.arcsec)
        b[refidx] = bfluxes[bidx]
        
        yra1 = self.yframe.stars['alphawin_j2000']
        ydec1 = self.yframe.stars['deltawin_j2000']
        yfluxes = self.yframe.stars['flux_auto']/self.yframe.expt
        ycoords = SkyCoord(ra=yra1*u.degree, dec=ydec1*u.degree)
        
        yidx, refidx, _, _ = search_around_sky(ycoords, refcoord, 0.6*u.arcsec)
        y[refidx] = yfluxes[yidx]
        

        oby = -2.5*np.log10(b / y)
          
        offset = np.nanmean(oby-cby)
        print offset
        plt.axhline(offset, linestyle='--', color='g')
        plt.scatter(cby, oby-cby, edgecolor='none', alpha=0.5)
        plt.grid()
        plt.show()

        plt.close()

        
if __name__ == '__main__':
    cs = CalStromgren()
    #ref = refstars()
    print cs.reference.stars['starid']
    cs.match()
    #print ref['starid']
    