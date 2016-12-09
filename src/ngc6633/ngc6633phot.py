#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 30, 2016

@author: Joerg Weingrill <jweingrill@aip.de>

reference stars:  1928
objids:  405
Cleaning ... 405 objids
Cleaning ... 1928 stars
(1928, 405)
(1928, 398)
(679, 398)
deleted  7 objids:  398 remaining
deleted  1249 stars:  679 remaining
(679, 398)
398 679
398 679
(557, 293) 557 557
(557, 293)

'''

import config
from datasource import DataSource
from cairo._cairo import Matrix
wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
import matplotlib.pyplot as plt
import numpy as np

def fprint(s):
    with open(config.datapath+'/ngc6633phot.log', 'at') as f:
        f.write(s)

def calibrate(objid, filtercol):
    ffield = {'V': '"V"', 'B': '"B"', 'I': '"I"'}
    queryparams = {'objid': objid, 
                   'ffield': ffield[filtercol]}
    query = """SELECT mag_isocor, magerr_isocor, %(ffield)s FROM phot, referencestars 
                WHERE objid = '%(objid)s'
                AND circle(phot.coord,0) <@ circle(referencestars.coord, 1.0/3600.0)
                AND flags < 4 AND mag_isocor>0 AND %(ffield)s > 0""" % queryparams
    stars = wifsip.query(query)
    if len(stars) < 15:
        print '0 0 0 0 0 0'
        fprint('0 0 0 0 0 0\n')
        return
    mag_isocor = np.array([r[0] for r in stars])
    #magerr_isocor = np.array([r[1] for r in stars])
    refmag = np.array([r[2] for r in stars])
    ocmag = mag_isocor - refmag
    std = np.std(ocmag)
    print '%d %f %f' % (len(stars), np.median(ocmag), std),
    fprint('%d %f %f ' % (len(stars), np.median(ocmag), std))
    mag_isocor1 = mag_isocor - np.median(ocmag)
    ocmag1 = mag_isocor1 - refmag
    k = abs(ocmag1 - std) < std
    if len(k) < 3:
        print '0 0 0'
        fprint('0 0 0\n')
        return
    mag_isocor2 = mag_isocor[k]
    refmag2 = refmag[k]
    mean =  np.mean(mag_isocor2 - refmag2)
    
    print '%d %f %f' % (len(mag_isocor2), mean, np.std(mag_isocor2 - refmag2))
    fprint('%d %f %f\n' % (len(mag_isocor2), mean, np.std(mag_isocor2 - refmag2)))
    #plt.title(objid)
    #plt.plot(refmag2, mag_isocor2-refmag2-mean, 'ob')
    #plt.xlabel('%s' % ffield[filtercol])
    #plt.show()

import pickle

class Calibrate(object):
    
    def __init__(self, filtercol):
        if filtercol not in ['B','V']:
            raise(ValueError,'Wrong filter color')
        self.filtercol= filtercol
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        
    def getreferencestars(self):
        query = 'SELECT starid, ra, dec, "B", "V" ' + \
                ' FROM referencestars ' + \
                ' WHERE referencestars.coord <@ box(point(276.44583,6.1416),point(277.17916,6.8749)) ' + \
                ' AND "V">0 and "B">0;'
        stars = self.wifsip.query(query)
        
        self.starids = [r[0] for r in stars]      
        self.ra = [r[1] for r in stars]      
        self.dec = [r[2] for r in stars] 
        if self.filtercol =='B':     
            self.refmag = np.array([r[3] for r in stars])      
        if self.filtercol =='V':     
            self.refmag = np.array([r[4] for r in stars])      
        self.safepickle('starids%s.pickle' % self.filtercol, self.starids)  
        print 'reference stars: ', self.nrefstars
    
    def getobjids(self):
        
        query = "SELECT objid FROM frames WHERE object LIKE 'NGC 6633 BVI %%' AND filter='%s'" % self.filtercol
        result = self.wifsip.query(query)
        self.objids = [r[0] for r in result]  
        self.safepickle('objids%s.pickle' % self.filtercol, self.objids)  
        print 'objids: ', self.nobjids
        self.mag = np.zeros((self.nrefstars, self.nobjids))
        self.err = np.zeros((self.nrefstars, self.nobjids))
    
    @property
    def nobjids(self):
        return len(self.objids)
    
    @property
    def nrefstars(self):
        return len(self.starids)
    
    def fillmatrix(self):
        try:
            self.mag = np.load(config.datapath+'mag%sarray.npy' % self.filtercol)
            self.err = np.load(config.datapath+'err%sarray.npy' % self.filtercol)
        except IOError:
            for starid, ra, dec in zip(self.starids, self.ra, self.dec):
                queryparams = {'filtercol': self.filtercol,
                               'ra': ra,
                               'dec': dec}
                query = """SELECT frames.objid, mag_isocor, magerr_isocor FROM phot, frames 
                    WHERE object LIKE 'NGC 6633 BVI %%' AND filter='%(filtercol)s'
                    AND phot.objid = frames.objid
                    AND circle(phot.coord,0) <@ circle(POINT(%(ra)f,%(dec)f), 1.0/3600.0)
                    AND flags < 4 AND mag_isocor>0""" % queryparams
                result = self.wifsip.query(query)
                objids = [r[0] for r in result]
                mags =  [r[1] for r in result]
                starindex = self.starids.index(starid)
                print starindex, starid, len(objids)
                for objid, mag, err in result:
                    objindex = self.objids.index(objid)
                    self.mag[starindex, objindex] = mag
                    self.err[starindex, objindex] = err
            np.save(config.datapath+'mag%sarray.npy' % self.filtercol, self.mag)
            np.save(config.datapath+'err%sarray.npy' % self.filtercol, self.err)
            self.array_toimage()

    def safepickle(self, filename, data):
        with open(config.datapath+filename, 'wb') as picklefile:
            pickle.dump(data, picklefile)

    def pickletonumpy(self, filename):
        with open(config.datapath+filename, 'rb') as picklefile:
            data = pickle.load(picklefile)
        newfilename = filename.replace('.pickle', '.npy')
        np.save(config.datapath+newfilename, data)
            
    def array_toimage(self, matrix=None, filename=None):
        from numpy import rint, isnan
        from PIL import Image  # @UnresolvedImport
        from functions import scaleto
        
        if matrix is None:
            m1 = self.mag
        else:
            m1 = matrix
            
        if filename is None:
            filename = 'mag%sarray.png'% self.filtercol
        
        m1[isnan(m1)] = 0.0
        simg = scaleto(m1,[0.0,255.0])
    
        simg = rint(simg)
        simg = simg.astype('uint8')
        im = Image.fromarray(simg)
        
        im.save(config.plotpath+filename) 
        
    def cleanmatrix(self, verbose=True):
        try:
            self.mag = np.load(config.datapath+'mag%sarray_clean_.npy' % self.filtercol)
            
        except IOError:
            delvec0 = []
            newobjids = []
    
            self.mag[self.mag<1.0] = np.nan
            
            print 'Cleaning ...',self.nobjids,'objids'
            
            
            #if an epoch shows less than 25% of the stars: remove it
            for i in range(self.nobjids):
                objidvec = self.mag[:, i]
                if verbose: print objidvec[np.isfinite(objidvec)].size,self.nrefstars
                if objidvec[np.isfinite(objidvec)].size < self.nrefstars/25:
                    delvec0.append(i)
                else:
                    newobjids.append(self.objids[i])
            #print delvec0
            
            
            print 'Cleaning ...',self.nrefstars,'stars'
            
            delvec1 = []
            newstarlist = []
            
            #if a star appears in less than 25% of the epochs: remove it 
            for j in range(self.nrefstars):
                starvec = self.mag[j, :]
                if verbose: 
                    print starvec[np.isfinite(starvec)].size, self.nobjids
                if starvec[np.isfinite(starvec)].size < self.nobjids/25:
                    delvec1.append(j)
                else:
                    newstarlist.append(list(self.starids)[j])
            
            
            print np.shape(self.mag)
            self.mag = np.delete(self.mag, delvec0, 1)
            print np.shape(self.mag)
            self.mag = np.delete(self.mag, delvec1, 0)
            print np.shape(self.mag)
            self.objids = newobjids
            self.refmag = np.delete(self.refmag, delvec1)
            print 'deleted ',len(delvec0),'objids: ', self.nobjids, 'remaining'
            self.starids = set(newstarlist)
            print 'deleted ',len(delvec1),'stars: ', self.nrefstars, 'remaining'
            print np.shape(self.mag)
            print self.nobjids, self.nrefstars
            assert np.shape(self.mag) == (self.nrefstars, self.nobjids)
            np.save(config.datapath+'mag%sarray_clean.npy' % self.filtercol, self.mag)                   
        
    def correctmatrix(self):
        
        inan = np.isnan(self.mag)
        refmatrix = np.tile(self.refmag, (self.nobjids, 1)).T
        assert(np.shape(refmatrix)==np.shape(self.mag))
        self.mag = self.mag - refmatrix
        self.array_toimage(filename='mag%sarray_corrected.png' % self.filtercol)
        #np.save(config.datapath+'mag%sarray_clean.npy' % self.filtercol, self.mag) 
        
        # restoring nan values that where broken after substraction
        self.mag[inan] = np.nan
        
        
        
        # remove objid outliers
        objidmean = np.nanmean(self.mag, axis=0)
        objidstd = np.nanstd(self.mag, axis=0)
        
        i = np.abs(objidmean - np.mean(objidmean))<np.std(objidmean)
        j = np.abs(objidstd - np.mean(objidstd))<np.std(objidstd)
        
        # remove refstars outliers
        
        starsmean = np.nanmean(self.mag, axis=1)
        starsstd = np.nanstd(self.mag, axis=1)
        
        k = np.abs(starsmean - np.mean(starsmean))<np.std(starsmean)
        l = np.abs(starsstd - np.mean(starsstd))<np.std(starsstd)
        print len(i & j), len(k & l)
        self.mag = self.mag[k & l, :]
        self.mag = self.mag[:, i & j]
        
        
        self.refmag = self.refmag[k & l]
        self.starids = np.array(list(self.starids))[k & l]
        self.objids = np.array(self.objids)[i & j]
        
        objidcorr = np.nanmean(self.mag, axis=0)
        print np.shape(self.mag), len(self.refmag), self.nrefstars
        corrmatrix = np.tile(objidcorr, (self.nrefstars, 1))
        print np.shape(corrmatrix)
        assert(np.shape(corrmatrix)==np.shape(self.mag))
        self.mag = self.mag - corrmatrix
        
        # remove refstars outliers
        
        starsstd = np.nanstd(self.mag, axis=1)
        
        l = starsstd<0.02
        self.mag = self.mag[l, :]
        self.refmag = self.refmag[l]
        self.starids = self.starids[l]
        
        # remove single outliers
        
        #magmean = np.nanmean(self.mag)
        #magstd = np.nanstd(self.mag)
        #m = np.abs(self.mag - magmean)> 3*magstd
        #self.mag[m] = np.nan
        
        plt.imshow(self.mag, interpolation='none', aspect='equal', vmin=-0.1, vmax=0.1)
        plt.show()
        plt.plot(self.refmag,np.nanstd(self.mag,axis=1),'or')
        plt.show()
        
        return
        
        

if __name__ == '__main__':
    cal = Calibrate('V')
    cal.getreferencestars()
    cal.getobjids()
    cal.fillmatrix()
    cal.cleanmatrix(verbose=False)
    cal.correctmatrix()
    exit()
    
    fprint('#objid filter expt field stars median std stars2 mean std2\n')
    for filtercol in ['V', 'B']:
        for expt in [60, 300]:
            if filtercol == 'B': expt *=2
            for field in ['C', 'NW', 'NE', 'SW', 'SE']:
                queryparams = {'filtercol': filtercol, 
                               'expt': str(expt),
                               'field': field}
                query = "SELECT objid FROM frames " + \
                        " WHERE object LIKE 'NGC 6633 BVI %(field)s' AND filter='%(filtercol)s' AND expt=%(expt)s;" % queryparams
                #print query
                objids = wifsip.query(query)
                for objid in objids:
                    print objid[0], filtercol, expt, field,
                    fprint('%s %s %d %s ' % (objid[0], filtercol, expt, field))
                    calibrate(objid[0], filtercol)
                
                