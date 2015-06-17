#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 2, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

perform differential photometry on the M48 dataset using PCA
'''
#import config
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

from diffphot import DiffPhotometry
class DiffPhotometryNoPCA(DiffPhotometry):
    '''
    inherited class from DiffPhotometry with different processing in reduced 
    and clean procedures
    '''
        
    def reduce(self):
        M = self._loadmatrix('phot')
        
        epochs, numstars = M.shape
        print 'reducing (%d,%d)' % M.shape
        print 'hjds: %d' % len(self.hjds)
        print 'objids: %d' % len(self.objids)
        print 'stars: %d' % len(self.stars)
        print 'starids: %d' % len(self.starids)
        
        
        
        assert len(self.stars) == len(self.starids)
        assert epochs == len(self.hjds)
        assert len(self.hjds) == len(self.objids)
        assert numstars == len(self.stars)
        
        maxnan = 1
        print 'removing stars and epochs'
        while maxnan>0.5:
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
                print 'e%3d, ' % delvec0,
            elif maxnan1 > maxnan0:
                M = np.delete(M, delvec1, 1)
                self.stars.pop(delvec1)
                self.starids.pop(delvec1)
                print 's%4d, ' % delvec1,
            else:
                print 'done (%.4f, %.4f)' % (maxnan0, maxnan1)
                
            maxnan = max(maxnan0,maxnan1)
            epochs, numstars = M.shape
        
        self._savematrix(M,'reduced')
        
        print 'saved reduced matrix', M.shape
        
        self._saveimage(self.filename+'_reducedmatrix.png', M, sort='mean')

    def clean(self, twosigma=False, sigmaclip=5.0):
        
        M = self._loadmatrix('reduced')
        
        print 'matrix to clean:', M.shape

        # calculate the mean for each lightcurve-vector
        meanvec = np.nanmean(M, axis=0)
        
        O = np.ones(M.shape)
        O = O * meanvec
        
        # calculate the mean for each epoch
        meanvec1 = np.nanmean(M-O, axis=1)
        
        P = np.ones(M.shape)
        P = (P.T * meanvec1).T
        
        M1 = M - O - P
        if twosigma:
            # create a new matrix with mean zero in time and in lightcurve
            # remove the stars where the stddev is larger than 2 sigmas
#            m1std = np.std(M1)*2.0
#            stdvec = np.where(np.std(M1, axis=0) > m1std)[0]
            
#            M1 = np.delete(M1, stdvec, 1)
#            O = np.delete(O, stdvec, 1)
#            P = np.delete(P, stdvec, 1)
#            if len(stdvec) > 1:
#                for i in stdvec: self.starids.pop(i)
#            elif len(stdvec) == 1: 
#                self.starids.pop(stdvec)
                
#            print 'matrix after 2sigma star reduction:', M1.shape
            
            m1std = np.nanstd(M1)*2.0
            stdvec = np.where(np.nanstd(M1, axis=1) > m1std)[0]
            
            fig, [ax1, ax2] = plt.subplots(1, 2)
            ax1.plot(np.sort(np.nanstd(M1, axis=0)))
            ax1.axhline(m1std)
            ax1.set_ylabel('std')
            ax1.set_xlabel('stars')
            
            ax2.plot(np.sort(np.nanstd(M1, axis=1)))
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
            else:
                print stdvec
                    
            print 'matrix after 2sigma epoch reduction:', M1.shape
        M = M1
        refstar = np.argmin(np.nanstd(M, axis=0))
        print 'reference star: %d, std = %.4f' % (refstar, np.nanstd(M[:, refstar]))
        
        # calculate where the std for each lightcurve-vector is < 0.01 mag
        stdvec = np.nanstd(M, axis=0)
        i = np.where(stdvec < 0.02)[0]
        print 'stars with std<0.02: %d' % len(i)
        # calculate the mean for each epoch
        meanvec2 = np.nanmean(M[:, i], axis=1)
        
        Q = np.ones(M.shape)
        Q =   (Q.T * meanvec2).T
        
        M = M - Q
        refstar = np.argmin(np.nanstd(M, axis=0))
        print 'new reference star: %d, std = %.4f' % (refstar, np.std(M[:, refstar]))

        # test for trends in matrix: ##########################################
        
        epochs, numstars = M.shape
        k = np.argsort(meanvec)
        w = 0.006*1.6**(meanvec[k]-14.3)
        R = M[:, k]
        for i in range(epochs):
            #plt.scatter(np.arange(numstars), R[i, :], edgecolor='none')
            x = np.arange(numstars)
            y = R[i, :]
            j = np.isfinite(y)
            z = np.polyfit(x[j], y[j], 2, w=1./w[j])
            if i==9:
                plt.scatter(x[j], y[j], edgecolor='none',alpha=0.5)
                plt.plot(x, np.polyval(z, x),'r')
                plt.show()
            R[i, :] -= np.polyval(z, np.arange(numstars))
            
        plt.imshow(R)
        plt.show()
        self._saveimage(self.filename+'_cleanedmatrixR.png', R, sort='std')
        
        plt.semilogy(meanvec[k], np.nanstd(R, axis=0), 'ro')
        plt.semilogy(meanvec,stdvec, 'o', mec='g', fillstyle='none')
        plt.xlabel('V mag')
        plt.ylabel('sigma')
        plt.grid()
        plt.savefig('/work2/jwe/owncloud/M48/plots/cleanedmatrixR.pdf')
        #plt.show()
        
        
        # ######################################################################
        
        self._savematrix(M, 'cleaned')
        self._saveimage(self.filename+'_cleanedmatrix.png', M, sort='std')
        
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
            filename = self.lightcurvepath+'/%s.dat' % starid
            a[:, 0] = self.hjds
            a[:, 1] = M[:, i]
            np.savetxt(filename, a, fmt='%.5f %.4f')
        
    def _plot_lightcurve(self, starid, t, m, axis):
        """
        plot the lightcurve for a given star
        """
        m -= np.mean(m)
        t -= min(t)
        plt.xticks(np.arange(0,100,10))
        plt.yticks(np.arange(-0.5,0.5,0.01))
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
                    plt.savefig(self.plotpath+self.field+'_lc%d.pdf' % lc)
                lc += 1
                sp = 1 
                plt.close()

