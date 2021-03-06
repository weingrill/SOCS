#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Mar 7, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

from frame import Frame
import numpy as np
from matplotlib import pyplot as plt
#from astropy.units import cbyte

def refstars():
    from datasource import DataSource
    
    table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
    query = "SELECT * FROM referencestars order by starid;"
    data = table.query(query)
    
    columns = data[0].keys()
    dtypes = [type(d) for d in data[0].values()]
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
            dtypes.append(type(d0[c]))
        
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
        '''
        20151126A-0008-0007 | M 67 uvby C |  480 | v
        20151126A-0008-0008 | M 67 uvby C |  600 | u
       
        '''
        
        self.bframe = Frame('20151126A-0008-0006')    
        self.yframe = Frame('20151126A-0008-0005') 
        self.vframe = Frame('20151126A-0008-0007')
        self.uframe = Frame('20151126A-0008-0008')
        
        self.hbnframe = Frame('20140227A-0011-0006') #'20151127A-0011-0004'
        self.hbwframe = Frame('20140227A-0011-0003') #'20151127A-0011-0003'
           
        
    def loadref(self):
        from astropy.coordinates import SkyCoord   # @UnresolvedImport
        from astropy import units as u

        self.reference = ReferenceStars()

        self.nstars = len(self.reference.stars['ra'])
        ref_ra = self.reference.stars['ra']
        ref_dec = self.reference.stars['dec']
        self.refcoord = SkyCoord(ra=ref_ra*u.degree, dec=ref_dec*u.degree)
    
    def refmatch(self, sourceframe):
        from astropy.coordinates import SkyCoord,search_around_sky   # @UnresolvedImport
        from astropy import units as u

        ra1 = sourceframe.stars['alphawin_j2000']
        dec1 = sourceframe.stars['deltawin_j2000']
        coords = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
        fluxes = sourceframe.stars['flux_auto']/sourceframe.expt
        
        bidx, refidx, _, _ = search_around_sky(coords, self.refcoord, 0.6*u.arcsec)
        # create an empty arry and set with nans
        compfluxes = np.ones(self.nstars) * np.nan
        
        compfluxes[refidx] = fluxes[bidx]
        
        return compfluxes
    
    def calbeta(self):
        cbeta = self.reference.stars['beta']
        cby = self.reference.stars['b-y']
        
        hbn = self.refmatch(self.hbnframe)
        hbw = self.refmatch(self.hbwframe)
        
        # convert to magnitudes
        obeta = -2.5*np.log10(hbn/hbw)
        
        # remove NaN values
        i = np.where(np.isfinite(obeta-cbeta))[0]

        cbeta = cbeta[i]
        obeta = obeta[i]
        cby   = cby[i]
        
        ones = np.ones(len(cbeta))
        
        A = np.vstack([ones, obeta, cby]).T
        c, kbeta, kby = np.linalg.lstsq(A, cbeta)[0]
        print('c = %.4f; kbeta = %.4f; kby = %.4f' %(c, kbeta, kby))
        
        o1beta = kbeta * obeta + kby * cby + c
        #betaoffset = np.nanmean(obeta-cbeta)
        
        #plt.axhline(betaoffset, linestyle='--', color='g')
        plt.plot(cby, obeta, 'bo', alpha=0.5)
        plt.plot(cby, o1beta, 'ro', alpha=0.5)
        plt.plot(cby, cbeta, 'go', alpha=0.5)
        
        plt.xlabel('beta')
        plt.ylabel('o - c')
        plt.grid()
        
        plt.show()
    
    def leastsqcal(self):
        cby = self.reference.stars['b-y']
        m1 = self.reference.stars['m1']
        c1 = self.reference.stars['c1']
        
        u = self.refmatch(self.uframe)
        v = self.refmatch(self.vframe)
        b = self.refmatch(self.bframe)
        y = self.refmatch(self.yframe)
        
        i = np.where(np.isfinite(u+v+b+y+cby+m1+c1))[0]
        
        cby = cby[i]
        m1 = m1[i]
        c1 = c1[i]
        u = -2.5*np.log10(u[i])
        v = -2.5*np.log10(v[i])
        b = -2.5*np.log10(b[i])
        y = -2.5*np.log10(y[i])
        
        ones = np.ones(len(cby))
        
        A = np.vstack([ones, b, y]).T
        c, kb, ky = np.linalg.lstsq(A, cby)[0]
        print(kb, ky, c)

        oby = kb*b + ky*y +c
        
        A = np.vstack([ones, v, b, y]).T
        c, kv, kb, ky = np.linalg.lstsq(A, m1)[0]
        print('m1', kv, kb, ky, c)
        om1 = kv*v + kb*b + ky*y + c
        
        c, ku, kv, kb = np.linalg.lstsq(A, c1)[0]
        print('c1', ku, kv, kb, c)
        oc1 = ku*u + kv*v + kb*b + c
        
        plt.plot(cby, oby, 'bo')
        plt.plot(m1, om1, 'mo')
        plt.plot(c1, oc1, 'co')
        plt.plot([0.0,1.0],[0.0,1.0], linestyle='--')
        plt.grid()
        plt.show()
    
    def match(self):
        
        # fill in the coordinates and b-y
        cby = self.reference.stars['b-y']
        cvb = self.reference.stars['m1'] + cby
        cuv = self.reference.stars['c1'] + cvb
        cbeta = self.reference.stars['beta']
        
        u = self.refmatch(self.uframe)
        v = self.refmatch(self.vframe)
        b = self.refmatch(self.bframe)
        y = self.refmatch(self.yframe)
        hbn = self.refmatch(self.hbnframe)
        hbw = self.refmatch(self.hbwframe)
        
        oby = -2.5*np.log10(b / y)
        ouv = -2.5*np.log10(u / v)
        ovb = -2.5*np.log10(v / b)
        obeta = -2.5*np.log10(hbn / hbw)
          
        byoffset = np.nanmean(oby-cby)
        uvoffset = np.nanmean(ouv-cuv)
        vboffset = np.nanmean(ovb-cvb)
        betaoffset = np.nanmean(obeta-cbeta)

        print('b - y: %.4f %.4f' % (byoffset, np.nanstd(oby - cby)))
        print('u - v: %.4f %.4f' % (uvoffset, np.nanstd(ouv - cuv)))
        print('v - b: %.4f %.4f' % (vboffset, np.nanstd(ovb - cvb)))
        print('beta : %.4f %.4f' % (betaoffset, np.nanstd(obeta - cbeta)))

        fig = plt.figure()

        ax1 = fig.add_subplot(411)

        plt.axhline(byoffset, linestyle='--', color='g')
        ax1.scatter(cby, oby-cby, edgecolor='none', alpha=0.5)
        ax1.set_xlabel('b - y')
        ax1.set_ylabel('o - c')
        plt.grid()
        
        ax2 = fig.add_subplot(412)
        plt.axhline(uvoffset, linestyle='--', color='g')
        ax2.scatter(cuv, ouv-cuv, edgecolor='none', alpha=0.5)
        ax2.set_xlabel('u - v')
        ax2.set_ylabel('o - c')
        plt.grid()

        ax3 = fig.add_subplot(413)
        plt.axhline(vboffset, linestyle='--', color='g')
        ax3.scatter(cvb, ovb-cvb, edgecolor='none', alpha=0.5)
        ax3.set_xlabel('v - b')
        ax3.set_ylabel('o - c')
        plt.grid()

        ax4 = fig.add_subplot(414)
        plt.axhline(betaoffset, linestyle='--', color='g')
        ax4.scatter(cbeta, obeta, edgecolor='none', alpha=0.5)
        ax4.set_xlabel('beta')
        ax4.set_ylabel('o - c')
        plt.grid()
        
        plt.show()

        plt.close()

        
if __name__ == '__main__':
    cs = CalStromgren()
    #cs.match()
    cs.calbeta()
    #cs.leastsqcal()
    
    