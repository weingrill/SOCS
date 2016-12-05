#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 11, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

We stardardized our photometry on the Stetson system (Stetson 2000). We 
discovered that the astrometry between our stars and Stetson's were off by 
0.3 arcseconds in both directions. By comparison with UCAC4, we found that our 
astrometric positions are within 0.2 arcsec of the ones given by UCAC4, wheras 
the positions from Stetson showed a shift. Since Stetson provided the plate 
coordinates in x and y, a recalculation of the plate solution was possible.

'''
import config
import numpy as np

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
    
    def readpos(self, filename='NGC6633.pos'):
        posfile = open(filename, 'rt')
        data = np.genfromtxt(posfile, dtype=None, names=('RA','DEC','HH','MM', 'SS','DD','DM','DS','dX','dY', 'X', 'Y', 'ID'))
        for d in data:
            #print ">%(ID)s<" % d
            query = """UPDATE ngc6633ref 
            SET ra=%(RA).11f,dec=%(DEC).11f, coord=point(%(RA).11f,%(DEC).11f)  
            WHERE starid='%(ID)s';""" % d
            print query
            self.wifsip.execute(query)
        posfile.close()
    
    def setucac4(self, filename = 'Stetson_UCAC4.fit'):
        import pyfits
        
        hdulist = pyfits.open(filename)
        ucac = hdulist[1].data
        hdulist.close()
        print ucac.columns
        print ucac['ucacid']
        
        #clear existing coordinates
        self.wifsip.execute('UPDATE ngc6633ref SET coord=NULL;')
        
        ra = ucac['raj2000']
        dec = ucac['dej2000']
        x = ucac['X']
        y = ucac['Y']
        A = np.vstack([x, y, np.ones(len(x))]).T
        wcsra = np.linalg.lstsq(A, ra)[0]
        wcsdec = np.linalg.lstsq(A, dec)[0]
        print wcsra
        print wcsdec
        dx1, dy1, cx = wcsra
        dx2, dy2, cy = wcsdec
        param= {'dx1': dx1,
                'dy1': dy1,
                'cx': cx,
                'dx2': dx2,
                'dy2': dy2,
                'cy': cy}
        query = """UPDATE ngc6633ref SET coord=point(X*%(dx1).11g+Y*%(dy1).11g+%(cx).11f, X*%(dx2).11g+Y*%(dy2).11g+%(cy).11f);""" % param
        print query
        self.wifsip.execute(query)

       
    def calibrate(self):
        from tools import log
        ffield = {'V': 'vmag', 'B': 'bmag', 'I': 'imag'}
        for filtercol in ['V','B','I']:
            query = """SELECT objid 
            FROM frames 
            WHERE object LIKE 'NGC 6633 BVI %%' 
            AND filter = '%s';""" % filtercol
            frames = self.wifsip.query(query)
            for frame in frames:
                params = {'objid': frame[0],
                          'filterfield': ffield[filtercol]}
                query = """SELECT starid, mag_isocor, magerr_isocor, %(filterfield)s, mag_isocor - %(filterfield)s
                FROM phot, ngc6633ref
                WHERE objid = '%(objid)s'
                AND circle(phot.coord,0) <@ circle(ngc6633ref.coord, 1.0/3600.0)
                AND NOT %(filterfield)s IS NULL
                AND flags<4;""" % params
                result = self.wifsip.query(query)
                if len(result)>0:
                    ocs = np.array([r[4] for r in result])
                    errs = np.array([r[2] for r in result])
                    weights=1/errs
                    std = np.std(ocs)
                    i = np.where(abs(ocs-np.mean(ocs))<std)[0]
                    #for r in result:
                    #    log(config.logfile, '%-11s %.4f %.4f %.3f %6.3f' % r)
                    corr = 0.0
                    if len(i) > 4:
                        corr = np.average(ocs[i], weights=weights[i])
                        std = np.std(ocs[i])
                    log(config.logfile, '%-19s %1s %-7.4f %.3f %3d %3d' % \
                        (frame[0], filtercol, corr, std, len(i), len(ocs)))
                    params['corr'] = corr
                    params['good'] = 'NULL'
                    if std < 0.05:
                        params['good'] = 'TRUE'
                    else:
                        params['good'] = 'FALSE'
                    query = """UPDATE frames SET corr=%(corr)f, good=%(good)s WHERE objid ='%(objid)s';""" % params
                else:
                    query = """UPDATE frames SET corr=NULL, good=FALSE WHERE objid ='%(objid)s';""" % params
                self.wifsip.execute(query)

    def calibratebv(self):
        '''drop view ngc6633match;
        create view ngc6633match as 
            SELECT ngc6633ref.starid, ngc6633.vmag, ngc6633.bmag, ngc6633ref.vmag "vmag_ref", ngc6633ref.bmag "bmag_ref" 
            FROM ngc6633, ngc6633ref 
            WHERE circle(ngc6633.coord,0) <@ circle(ngc6633ref.coord, 1.0/3600.0);'''
        
        from dbtable import DBTable
        match = DBTable(self.wifsip, 'ngc6633match', condition = 'NOT vmag_ref IS NULL AND NOT bmag_ref IS NULL AND vmag>12')
        vmag = match['vmag']
        bmag = match['bmag']
        bv= bmag - vmag
        vmagref = match['vmag_ref']
        bmagref = match['bmag_ref']
        bvref = bmagref - vmagref
        A = np.vstack([vmag, bv, np.ones(len(vmag))]).T
        vcorr = np.linalg.lstsq(A, vmagref)[0]
        B = np.vstack([bmag, bv, np.ones(len(bmag))]).T
        bcorr = np.linalg.lstsq(B, bmagref)[0]
        
        print vcorr
        print bcorr
        

stet = Stetson()
#stet.readpos(filename=config.datapath+'NGC6633.pos')
#stet.setucac4(filename=config.datapath+'Stetson_UCAC4.fit')
#stet.fromfile('/work2/jwe/SOCS/NGC6633/data/NGC6633.fit')
#stet.todatabase()        
#stet.calibrate()
stet.calibratebv()