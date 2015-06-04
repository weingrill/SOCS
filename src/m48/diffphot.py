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

class Epoch(object):
    def __init__(self):
        self.frame_objid = ''
        self.hjd = 0.0
    
class Star(object):
    def __init__(self, starid):
        self.starid = starid
        self.mag = 0.0
        self.coord = (0.0, 0.0)
    
    @property    
    def objid(self):
        return self.starid.split('#')[0]

    @property    
    def star(self):
        return self.starid.split('#')[1]


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
        #self.starids = []
        self.objids = []
        self.field = field
        self.filename = config.datapath+self.field
        print self.field
        
    
    def load_objids(self):
        """
        load the objids corresponding to the field
        """
        query = """SELECT objid, hjd, stars 
        FROM frames 
        WHERE object LIKE '%s' 
        AND filter='V'
        AND expt>60
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
        print 'maximum stars: %4d' % self.numstars
    
    def _loadcoords(self, objid):
        query = """SELECT star, coord 
        FROM phot 
        WHERE objid = '%s'
        AND flags<4
        ORDER BY star""" % objid
        result = self.wifsip.query(query)
        print 'number of coordinates: %4d' % len(result)
        stars = [r[0] for r in result]
        coords = [r[1] for r in result] 
        return stars, coords
    
    def _loadmags(self, coords):
        #print '%s'  % coords
        query = """SELECT frames.objid, mag_auto 
            FROM phot, frames
            WHERE object like '%s'
            AND frames.objid = phot.objid
            AND filter='V'
            AND expt > 60
            AND phot.flags<4
            AND circle(phot.coord,0) <@ circle(point%s, 0.2/3600.0)""" % \
            (self.field, coords)
        result = self.wifsip.query(query)
        objids = [r[0] for r in result] 
        mags = np.array([r[1] for r in result])
        return objids, mags
    
    def _savestars(self):
        f = open(self.filename+'_starids.txt', 'wt')
        
#         for s in self.stars:
#             f.write('%s#%d\n' % (self.refobjid, s))
        
        for s in self.starids:
            f.write('%s\n' % s)
        f.close()

        
    def _saveepochs(self):
        f = open(self.filename+'_epochs.txt', 'wt')
        for o, h in zip(self.objids, self.hjds):
            f.write('%s %f\n' % (o, h))
        f.close()
        
    def build_photmatrix(self):
        """
        creates the photmatrix by loading all stars for each objid/epoch
        """
        print 'reference objid: %s' % self.refobjid
        
        self.stars, refcoords = self._loadcoords(self.refobjid)
        self.numstars = len(self.stars)
        self.starids = ['%s#%d' % (self.refobjid, s) for s in self.stars]
        self._savestars()

        print 'number of reference stars: %d' % self.numstars
        
        try:
            photmatrix = np.load(self.filename+'_photmatrix.npy')
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
        
        epochs, numstars = photmatrix.shape
        assert epochs == len(self.hjds)
        assert numstars == len(self.stars)
        np.save(self.filename+'_photmatrix.npy', photmatrix)
        self._savestars()
        self._saveepochs()
        self._saveimage(self.filename+'_photmatrix.png', photmatrix)
        
    
    def _saveimage(self, filename, myarray):
        """
        
        First ensure your numpy array, myarray, is normalised with the max value at 1.0.
        Apply the colormap directly to myarray.
        Rescale to the 0-255 range.
        Convert to integers, using np.uint8().
        Use Image.fromarray().

        """
        from PIL import Image  # @UnresolvedImport
        
        myarray *= -1.0
        myarray -= np.nanmin(myarray)
        myarray *= np.trunc(256.0/np.nanmax(myarray))
        
        simg = np.rint(myarray)
        simg = simg.astype('uint8')
        
        im = Image.fromarray(plt.cm.jet(simg, bytes=True))  # @UndefinedVariable
        im.save(filename)
        
        #plt.imshow(myarray, interpolation='nearest', cmap=cmap)
        #plt.axis('off')
        #plt.savefig(filename)
        #plt.close()
        
    
    def reduce(self):
        try:
            M = np.load(self.filename+'_reducedmatrix.npy')
        except IOError:
            M = np.load(self.filename+'_photmatrix.npy')
            #self.load_objids()
            #self.stars, _ = self._loadcoords(self.refobjid)
        else:
            return
        
        try:
            self._loadstars()
        except IOError:
            self.stars, _ = self._loadcoords(self.refobjid)
            self.starids = ['%s#%d' % (self.refobjid, s) for s in self.stars]
        try:
            self._loadepochs()
        except IOError:
            self.load_objids()
        #
        
        epochs, numstars = M.shape
        print 'reducing (%d,%d)' % M.shape
        print 'hjds, objids, stars = %d, %d' % (len(self.hjds), len(self.stars))
        
        assert len(self.stars) == len(self.starids)
        assert epochs == len(self.hjds)
        assert len(self.hjds) == len(self.objids)
        assert numstars == len(self.stars)
        
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
                self.objids.pop(delvec0)
                self.hjds = np.delete(self.hjds, delvec0)
                print 'deleted epoch %3d (%.3f)' % (delvec0, maxnan0)
            elif maxnan1 > maxnan0:
                M = np.delete(M, delvec1, 1)
                self.stars.pop(delvec1)
                self.starids.pop(delvec1)
                print 'deleted star %4d (%.3f)' % (delvec1, maxnan1)
            else:
                print maxnan0, maxnan1
                
            maxnan = max(maxnan0,maxnan1)
            epochs, numstars = M.shape
        
        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)
        np.save(config.datapath+field+'_reducedmatrix.npy', M)
        self._saveepochs()
        self._savestars()
        
        print 'saving reduced matrix', M.shape
        
        self._saveimage(self.filename+'_reducedmatrix.png', M)
        
    def clean(self, twosigma=False, sigmaclip=5.0):
        from scipy.linalg import svd
        
        try:
            M = np.load(self.filename+'_cleanedmatrix.npy')
        except IOError:
            M = np.load(self.filename+'_reducedmatrix.npy')
        else:
            return
        
        print 'matrix to clean:', M.shape
        self._loadepochs()
        self._loadstars()

        epochs, numstars = M.shape
        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)
        
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
        M1 = M - O #- P
        
#         plt.subplot(3, 2, 1)
#         plt.imshow(M)
#         plt.subplot(3, 2, 3)
#         plt.imshow(O)
#         plt.subplot(3, 2, 5)
#         plt.imshow(P)
#         plt.subplot(3, 2, 4)
#         v = 3.0
#         plt.imshow(M-O, vmin = -v, vmax = v)
#         plt.subplot(3, 2, 6)
#         plt.imshow(M-O-P, vmin = -v, vmax = v)
#         plt.show()
        
        if twosigma:
            # remove the stars where the stddev is larger than 2 sigmas
            m1std = np.std(M1)*2.0
            stdvec = np.where(np.std(M1, axis=0) > m1std)[0]
            
            M1 = np.delete(M1, stdvec, 1)
            O = np.delete(O, stdvec, 1)
            P = np.delete(P, stdvec, 1)
            if len(stdvec) > 1:
                for i in stdvec: self.starids.pop(i)
            elif len(stdvec) == 1: 
                self.starids.pop(stdvec)
                
            print 'matrix after 2sigma star reduction:', M1.shape
            
            m1std = np.std(M1)*2.0
            stdvec = np.where(np.std(M1, axis=1) > m1std)[0]
            
            fig, [ax1, ax2] = plt.subplots(1, 2)
            ax1.plot(np.sort(np.std(M1, axis=0)))
            ax1.axhline(m1std)
            ax1.set_ylabel('std')
            ax1.set_xlabel('stars')
            
            ax2.plot(np.sort(np.std(M1, axis=1)))
            ax2.axhline(m1std)
            ax2.set_xlabel('epochs')
            fig.savefig(self.filename+'_matrixstat.pdf')
            plt.close()
            
            # remove the epochs where the stddev is larger than 2 sigmas
            M1 = np.delete(M1, stdvec, 0)
            O = np.delete(O, stdvec, 0)
            P = np.delete(P, stdvec, 0)
            self.hjds = np.delete(self.hjds, stdvec)
            if len(stdvec) > 1:
                for i in stdvec: 
                    self.objids.pop(i)
            elif len(stdvec) == 1: 
                    self.objids.pop(stdvec)
                    
            print 'matrix after 2sigma epoch reduction:', M1.shape
            
        print 'reference star: ',np.argmin(np.std(M1, axis=0))
        
        
        # perform sigma clipping at a level of 3 sigmas
        m1std = np.std(M1)
        
        i = np.where(M1 > sigmaclip*m1std)
        M1[i] = sigmaclip*m1std
        j = np.where(M1 < -sigmaclip*m1std)
        M1[j] = -sigmaclip*m1std
        
        M = M1 + O
        epochs, numstars = M.shape
        
#         plt.subplot(3, 1, 1)
#         plt.imshow(M1)
#         plt.subplot(3, 1, 2)
#         plt.imshow(O)
#         plt.subplot(3, 1, 3)
#         plt.imshow(P)
#         plt.show()
        
        
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

        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)
        np.save(self.filename+'_cleanedmatrix.npy', M-Mhat2)
        self._saveepochs()
        self._savestars()
        plt.semilogy(s)
        plt.title('singular values')
        plt.savefig(self.filename+'_cleanedmatrixs.pdf')
        plt.close()
        self._saveimage(self.filename+'_cleanedmatrix.png', M-Mhat2)
        self._saveimage(self.filename+'_cleanedmatrixM.png', M)
        self._saveimage(self.filename+'_cleanedmatrixM1.png', M1)
        
    def _loadstars(self):
        f = open(self.filename+'_starids.txt', 'rt')
        self.starids = [starid.rstrip('\n') for starid in f.readlines()]
        f.close()
        
    def _loadepochs(self):
        f = open(self.filename+'_epochs.txt', 'rt')
        lines = f.readlines()
        f.close()
        
        self.objids = [l.split()[0] for l in lines]
        self.hjds = np.array([float(l.split()[1].rstrip('\n')) for l in lines])
        
    def save_lightcurves(self):
        M = np.load(self.filename+'_cleanedmatrix.npy')
        epochs, numstars = M.shape
        self._loadstars()
        self._loadepochs()
        
        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)
        
        a = np.empty([epochs, 2])
        
        for starid in self.starids:
            i = self.starids.index(starid)
            filename = config.projectpath+'lightcurves.new/%s.dat' % starid
            a[:, 0] = self.hjds
            a[:, 1] = M[:, i]
            np.savetxt(filename, a, fmt='%.5f %.4f')
        
    def _plot_lightcurve(self, starid, t, m, axis):
        """
        plot the lightcurve for a given star
        """
        plt.xlim(min(t),max(t))
        plt.scatter(t, m, edgecolor='none', facecolor='g', s=5)
        plt.plot(t,m,'gray')
        std = 5.0*np.std(m)
        plt.text(0.95, 0.95, '#'+starid.split('#')[1], 
                 fontsize=12,
                 verticalalignment='top',
                 horizontalalignment='right',
                 transform=axis.transAxes)
        plt.ylim([std, -std])

    def make_lightcurves(self, show=False):
        """plot lightcurve"""
        from matplotlib import rcParams
        
        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)
        
        M = np.load(self.filename+'_cleanedmatrix.npy')
        epochs, numstars = M.shape
        self._loadstars()
        self._loadepochs()
        
        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)
        
        sp = 1
        lc = 1
        rows = 7
        cols = 4
        for starid in self.starids:
            i = self.starids.index(starid)
            ax = plt.subplot(rows, cols,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self._plot_lightcurve(starid, self.hjds, M[:, i], ax)
            if sp == rows*cols + 1 or starid == self.starids[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig(config.plotpath+'newlightcurves%d.pdf' % lc)
                lc += 1
                sp = 1 
                plt.close()
        

if __name__ == '__main__':
    
#    fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE','M 48 rot SW']
    fields =    ['M 48 rot NE']
    for field in fields:
        diffphot = DiffPhotometry(field)
        diffphot.load_objids()
        diffphot.build_photmatrix()
        diffphot.reduce()
        diffphot.clean(twosigma = False)
        #diffphot.save_lightcurves()
        diffphot.make_lightcurves(show=True)
