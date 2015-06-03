#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 2, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

perform differential photometry on the M48 dataset using PCA
'''
import config
import numpy as np
import matplotlib.pyplot as plt

class DiffPhotometry(object):
    '''
    classdocs
    '''

    def __init__(self, field):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.wifsip = DataSource(database='stella', 
                                 host='pera', 
                                 user='stella')
        self.coords = {}
        self.starids = []
        self.objids = []
        self.field = field
        print self.field
        
    
    def load_objids(self):
        """
        load the objids corresponding to the field
        """
        query = """SELECT objid, hjd, stars 
        FROM frames 
        WHERE object LIKE '%s' 
        AND filter='V'
        ORDER BY objid;""" % self.field
        result = self.wifsip.query(query)
        
        self.objids = [r[0] for r in result] 
        self.epochs = len(self.objids)
        
        self.hjds = np.array([r[1] for r in result])
        starcounts =  [r[2] for r in result]
        
        self.numstars = max(starcounts)
        i = np.argmax(starcounts)
        self.refobjid = self.objids[i]
        print 'number of objids: %3d' % len(self.objids)
        print 'number of stars: %4d' % self.numstars
    
    def _loadcoords(self, objid):
        query = """SELECT coord 
        FROM phot 
        WHERE objid = '%s'
        ORDER BY star""" % objid
        result = self.wifsip.query(query)
        print 'number of coordinates: %4d' % len(result)
        return [r[0] for r in result]
    
    def _loadmags(self, coords):
        #print '%s'  % coords
        query = """SELECT frames.objid, mag_auto 
            FROM phot, frames
            WHERE object like '%s'
            AND frames.objid = phot.objid
            AND filter='V'
            AND phot.flags<4
            AND circle(phot.coord,0) <@ circle(point%s, 0.2/3600.0)""" % \
            (self.field, coords)
        result = self.wifsip.query(query)
        objids = [r[0] for r in result] 
        mags = np.array([r[1] for r in result])
        return objids, mags
        
    def build_photmatrix(self):
        """
        creates the photmatrix by loading all stars for each objid/epoch
        """
        print 'reference objid: %s' % self.refobjid
        
        refcoords = self._loadcoords(self.refobjid)
        print 'number of reference coords: %d' % len(refcoords)
        
        try:
            photmatrix = np.load(config.datapath+self.field+'_photmatrix.npy')
        except IOError:
            photmatrix = np.zeros([self.epochs,self.numstars])*np.nan
        else:    
            return
        print 'generating photmatrix with dimensions %d %d' % photmatrix.shape
        
        for star in range(len(refcoords)):
            coord = refcoords[star]
            objids, mags = self._loadmags(coord)
            
            print '%4d %s %4d' % (star, coord, len(objids))
            
            for objid in objids:
                epoch = self.objids.index(objid)
                k = objids.index(objid)
                photmatrix[epoch, star] = mags[k]
        
        np.save(config.datapath+self.field+'_photmatrix.npy', photmatrix)
        
        cmap = plt.cm.get_cmap('jet')
        cmap.set_bad('grey')
        plt.imshow(photmatrix, interpolation='nearest', cmap=cmap)
        plt.show()
    
    def reduce(self):
        try:
            M = np.load(config.datapath+self.field+'_reducedmatrix.npy')
        except IOError:
            M = np.load(config.datapath+self.field+'_photmatrix.npy')
        else:
            return
        
        epochs, numstars = M.shape
        print 'reducing (%d,%d)', M.shape
        maxnan = 1
        while maxnan>0:
            maxnan0 = 0
            delvec0 = None
            for i in range(epochs):
                objidvec = M[i,:]
                nancount0 = 1.0 * objidvec[np.isnan(objidvec)].size / len(objidvec)
                if nancount0 > maxnan0:
                    maxnan0 = nancount0
                    delvec0 = i
                    
        
            maxnan1=0    
            delvec1 = None
            for j in range(numstars):
                starvec = M[:,j]
                nancount1 =  1.0 * starvec[np.isnan(starvec)].size / len(starvec)
                if nancount1 > maxnan1:
                    maxnan1 = nancount1
                    delvec1 = j
                    
            if maxnan0 >= maxnan1 and maxnan0 > 0:
                M = np.delete(M, delvec0, 0)
                print 'deleted epoch %3d (%.3f)' % (delvec0, maxnan0)
            elif maxnan1 > maxnan0:
                M = np.delete(M, delvec1, 1)
                print 'deleted star %4d (%.3f)' % (delvec1, maxnan1)
            else:
                print maxnan0, maxnan1
                
            maxnan = max(maxnan0,maxnan1)
            epochs, numstars = M.shape
        
        #np.save(config.datapath+field+'_reducedmatrix.npy', M)
        
        print M.shape
        cmap = plt.cm.get_cmap('jet')
        cmap.set_bad('grey')
        plt.imshow(M, interpolation='nearest', cmap=cmap)
        plt.show()
    
    def clean(self):
        from scipy.linalg import svd
        
        try:
            M = np.load(config.datapath+self.field+'_cleanedmatrix.npy')
        except IOError:
            M = np.load(config.datapath+self.field+'_reducedmatrix.npy')
        else:
            return

        # calculate the mean for each lightcurve-vector
        meanvec = np.mean(M, axis=0)
        
        O = np.ones(M.shape)
        O = O * meanvec
        #M = M - O
        
        # calculate the mean for each epoch
        meanvec1 = np.mean(M-O, axis=1)
        
        P = np.ones(M.shape)
        P =   (P.T * meanvec1).T
        
        # create a new matrix with mean zero in time and in lightcurve
        M1 = M - O - P
        
        m1std = np.std(M1)*2.0
        stdvec = np.where(np.std(M1, axis=0) > m1std)
        #plt.plot(np.sort(np.std(M1, axis=0)))
        #plt.axhline(m1std)
        #plt.show()
        
        # remove the stars where the stddev is larger than 2 sigmas
        M1 = np.delete(M1, stdvec, 1)
        O = np.delete(O, stdvec, 1)
        P = np.delete(P, stdvec, 1)
        print M1.shape
        
        m1std = np.std(M1)*2.0
        stdvec = np.where(np.std(M1, axis=1) > m1std)
        
        #plt.plot(np.sort(np.std(M1, axis=1)))
        #plt.axhline(m1std)
        #plt.show()
        print stdvec
        
        # remove the epochs where the stddev is larger than 2 sigmas
        M1 = np.delete(M1, stdvec, 0)
        O = np.delete(O, stdvec, 0)
        P = np.delete(P, stdvec, 0)
        #t = np.delete(t, stdvec)
        print M1.shape
        
        print 'reference star: ',np.argmin(np.std(M1, axis=0))
        
        
        m1std = np.std(M1)
        s = 3.0
        
        i = np.where(M1 > s*m1std)
        M1[i] = s*m1std
        j = np.where(M1 < -s*m1std)
        M1[j] = -s*m1std
        
        M = M1 + O
        epochs, numstars = M.shape
        
        U, s, Vt = svd(M, full_matrices=False)
        V = Vt.T
        
        # sort the PCs by descending order of the singular values (i.e. by the
        # proportion of total variance they explain)
        ind = np.argsort(s)[::-1]
        U = U[:, ind]
        s = s[ind]
        V = V[:, ind]
        
        # if we use all of the PCs we can reconstruct the noisy signal perfectly
        S = np.diag(s)
        Mhat = np.dot(U, np.dot(S, V.T))
        #print "Using all PCs, MSE = %.6G" %(np.mean((M - Mhat)**2))
        
        # if we use only the first n PCs the reconstruction is less accurate
        level = 1
        Mhat2 = np.dot(U[:, :level], np.dot(S[:level, :level], V[:,:level].T))
        #print "Using first %d PCs, MSE = %.6G" %(level, np.mean((M - Mhat2)**2))

        #np.save(config.datapath+field+'_cleanedmatrix.npy', M-Mhat2)
        
        cmap = plt.cm.get_cmap('jet')
        cmap.set_bad('grey')
        plt.imshow(M-Mhat2, interpolation='nearest', cmap=cmap)
        plt.show()

        exit()
        k = 450
        
        for k in range(numstars):
            fig, [ax1, ax2] = plt.subplots(2, 1)
            ax1.set_title('star %d' % k)
            t = np.arange(len(M1[:, k]))
            ax1.scatter(t, M1[:, k], edgecolor='none')
            ax1.axis('tight')
            
            
            ax2.scatter(t, M[:, k]-Mhat2[:, k], edgecolor='none')
            ax2.axis(ax1.axis())
            
            #ax3.semilogy(s[:10])
            lc0std = np.std(M1[:, k])
            lc1std = np.std(M[:, k]-Mhat2[:, k])
            ax2.set_xlabel('%.4f --> %.4f' % (lc0std, lc1std  ))
            plt.show()
        
        Mimg = M - Mhat2
        plt.imshow(Mimg, 
                   interpolation='nearest')
        plt.show()        

    def _lightcurve(self, starid):
        print 'lightcurve for star: %s' % starid
        coords = self.coords[starid]
        query = """SELECT frames.objid, hjd, mag_auto 
            FROM phot, frames
            WHERE object like '%s'
            AND frames.objid = phot.objid
            AND filter='V'
            AND phot.flags<4
            AND circle(phot.coord,0) <@ circle(point(%f, %f), 0.2/3600.0)""" % \
            (self.field, coords[0], coords[1])
        result = self.wifsip.query(query)
        objids = [r[0] for r in result] 
        hjds = np.array([r[1] for r in result]) 
        mags = np.array([r[2] for r in result])
        return objids, hjds, mags
    
        
if __name__ == '__main__':
    
#    fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE','M 48 rot SW']
    fields =    ['M 48 rot NE']
    for field in fields:
        diffphot = DiffPhotometry(field)
        diffphot.load_objids()
        diffphot.build_photmatrix()
        diffphot.reduce()
        diffphot.clean()
