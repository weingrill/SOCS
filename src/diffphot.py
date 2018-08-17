#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 2, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

perform differential photometry on the M48 dataset using PCA
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.linalg import svd


class Epoch(object):
    """
    class to represent an epoch in the dataset
    """

    def __init__(self):
        self.frame_objid = ''
        self.hjd = 0.0


class Star(object):
    """
    class to represent a star in the dataset
    """

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
    class for performing differential photometry on a rot dataset
    '''

    def __init__(self, field,
                 datapath,
                 lightcurvepath,
                 plotpath,
                 filtercol='V',
                 flaglimit=4,
                 starslimit=10000,
                 dettemp=-110.0):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.wifsip = DataSource(database='stella',
                                 host='pera',
                                 user='stella')
        self.coords = {}
        self.objids = []
        self.field = field
        self.filtercol = filtercol
        self.flaglimit = flaglimit
        self.starslimit = starslimit
        self.dettemp = dettemp
        self.filename = datapath + self.field
        self.lightcurvepath = lightcurvepath
        self.plotpath = plotpath
        print(self.field)
        print('filter:   %s' % self.filtercol)
        print('maxstars: %d' % self.starslimit)
        print('dettemp:  <%.1f' % self.dettemp)

    def load_objids(self):
        """
        load the objids corresponding to the field
        """

        data = {'object': self.field,
                'filtercol': self.filtercol,
                'stars': self.starslimit,
                'dettemp': self.dettemp}
        query = """SELECT objid, hjd, stars 
            FROM frames 
            WHERE object LIKE '%(object)s' 
            AND filter='%(filtercol)s'
            AND stars<%(stars)d
            AND dettemp<%(dettemp)f
            ORDER BY objid;""" % data
        # AND expt>60

        result = self.wifsip.query(query)

        self.objids = [r[0] for r in result]
        self.epochs = len(self.objids)

        self.hjds = np.array([r[1] for r in result])
        starcounts = [r[2] for r in result]

        self.numstars = max(starcounts)
        i = np.argmax(starcounts)
        self.refobjid = self.objids[i]
        print('number of objids: %3d' % len(self.objids))
        print('maximum stars:    %4d' % self.numstars)

    def _loadcoords(self, objid):
        """
        load the star number and its coordinate
        """
        query = """SELECT star, coord 
        FROM phot 
        WHERE objid = '%s'
        AND flags<%d
        ORDER BY star""" % (objid, self.flaglimit)
        result = self.wifsip.query(query)
        print('number of coordinates: %4d' % len(result))
        stars = [r[0] for r in result]
        coords = [r[1] for r in result]
        return stars, coords

    def _loadmags(self, coords):
        """
        returns the frame's objid and the magnitude for the stars at coordinate
        'coords' throughout all the frames
        """
        data = {'object': self.field,
                'filter': self.filtercol,
                'flaglimit': self.flaglimit,
                'stars': self.starslimit,
                'dettemp': self.dettemp,
                'coords': coords}
        query = """SELECT frames.objid, mag_auto 
            FROM phot, frames
            WHERE object like '%(object)s'
            AND frames.objid = phot.objid
            AND filter = '%(filter)s'
            AND phot.flags < %(flaglimit)d
            AND stars < %(stars)d
            AND dettemp < %(dettemp)f
            AND circle(phot.coord,0) <@ circle(point%(coords)s, 0.2/3600.0)""" % data

        # AND expt>60
        result = self.wifsip.query(query)
        objids = [r[0] for r in result]
        mags = np.array([r[1] for r in result])
        return objids, mags

    def _savestars(self, stage=''):
        """
        saves the starids in a textfile; requires the correct stage to be set
        """
        if not stage in ['phot', 'reduced', 'cleaned']:
            raise (TypeError)
        f = open(self.filename + '_%sstarids.txt' % stage, 'wt')
        for s in self.starids:
            f.write('%s\n' % s)
        f.close()

    def _saveepochs(self, stage=''):
        """
        saves the epochs (objids) and the corresponding hjd in a textfile; 
        requires the correct stage to be set
        """
        if not stage in ['phot', 'reduced', 'cleaned']:
            raise (TypeError)
        f = open(self.filename + '_%sepochs.txt' % stage, 'wt')
        for o, h in zip(self.objids, self.hjds):
            f.write('%s %f\n' % (o, h))
        f.close()

    def build_photmatrix(self):
        """
        creates the photmatrix by loading all stars for each objid/epoch
        """
        from tools import log

        print('reference objid: %s' % self.refobjid)

        self.stars, refcoords = self._loadcoords(self.refobjid)
        self.numstars = len(self.stars)
        self.starids = ['%s#%d' % (self.refobjid, s) for s in self.stars]

        print('number of reference stars: %d' % self.numstars)

        try:
            photmatrix = np.load(self.filename + '_photmatrix.npy')
        except IOError:
            photmatrix = np.zeros([self.epochs, self.numstars]) * np.nan
        else:
            return
        print('generating photmatrix with dimensions %d %d' % photmatrix.shape)

        for star in range(len(refcoords)):
            coord = refcoords[star]
            objids, mags = self._loadmags(coord)

            log(self.filename + 'loadmags.txt',
                logstring='%4d %s %4d' % (star, coord, len(objids)))
            for objid in objids:
                epoch = self.objids.index(objid)
                k = objids.index(objid)
                photmatrix[epoch, star] = mags[k]

        epochs, numstars = photmatrix.shape
        assert epochs == len(self.hjds)
        assert numstars == len(self.stars)
        self._savematrix(photmatrix, 'phot')
        self._saveimage(self.filename + '_photmatrix.png', photmatrix)

    def _loadmatrix(self, stage):
        """
        loads the matrix at the given stage and sets the following properties:
        objids
        hjds,
        stars,
        starids
        """
        if not stage in ['phot', 'reduced', 'cleaned']:
            raise (TypeError)
        matrix = np.load(self.filename + '_%smatrix.npy' % stage)

        try:
            f = open(self.filename + '_%sstarids.txt' % stage, 'rt')
            self.starids = [starid.rstrip('\n') for starid in f.readlines()]
            f.close()
        except IOError:
            self.load_objids()  # to populate self.refobjid
            self.stars, _ = self._loadcoords(self.refobjid)
            self.starids = ['%s#%d' % (self.refobjid, s) for s in self.stars]

        try:
            f = open(self.filename + '_%sepochs.txt' % stage, 'rt')
            lines = f.readlines()
            f.close()
            self.objids = [l.split()[0] for l in lines]
            self.hjds = np.array([float(l.split()[1].rstrip('\n')) for l in lines])
        except IOError:
            self.load_objids()

        self.stars = [int(starid.split('#')[1]) for starid in self.starids]
        epochs, numstars = matrix.shape
        assert len(self.stars) == len(self.starids)
        assert epochs == len(self.hjds)
        assert len(self.hjds) == len(self.objids)
        assert numstars == len(self.stars)

        return matrix

    def _savematrix(self, matrix, stage):
        """
        saves the matrix into a numpy file at the given stage
        """
        if not stage in ['phot', 'reduced', 'cleaned']:
            raise (TypeError)
        epochs, numstars = matrix.shape

        assert epochs == len(self.hjds)
        assert numstars == len(self.starids)

        np.save(self.filename + '_%smatrix.npy' % stage, matrix)

        self._savestars(stage)
        self._saveepochs(stage)

    def _saveimage(self, filename, myarray, sort='std'):
        """
        First ensure your numpy array, myarray, is normalised with the max value at 1.0.
        Apply the colormap directly to myarray.
        Rescale to the 0-255 range.
        Convert to integers, using np.uint8().
        Use Image.fromarray().
        """
        from PIL import Image  # @UnresolvedImport

        mavg = np.nanmean(myarray)

        myarray -= mavg
        mstd = 2.0 * np.nanstd(myarray)
        i = np.where(np.isnan(myarray))
        myarray[i] = 0.0
        if sort == 'std':
            s = np.nanstd(myarray, axis=0)
            ind = np.argsort(s)
            myarray = myarray[:, ind]
        elif sort == 'mean':
            s = np.nanmean(myarray, axis=0)
            ind = np.argsort(s)
            myarray = myarray[:, ind]

        i = np.where(myarray > mstd)
        myarray[i] = mstd
        j = np.where(myarray < -mstd)
        myarray[j] = -mstd

        myarray *= -1.0
        myarray -= np.nanmin(myarray)
        myarray *= np.trunc(256.0 / np.nanmax(myarray))

        simg = np.rint(myarray)
        simg = simg.astype('uint8')

        im = Image.fromarray(plt.cm.jet(simg, bytes=True))  # @UndefinedVariable
        im.save(filename)

        # plt.imshow(myarray, interpolation='nearest', cmap=cmap)
        # plt.axis('off')
        # plt.savefig(filename)
        # plt.close()

    def reduce(self):
        """
        reduce the original matrix by removing epochs or starswith the most nan 
        values successively
        """
        M = self._loadmatrix('phot')
        epochs, numstars = M.shape

        print('reducing (%d,%d)' % M.shape)

        maxnan = 1
        print('removing stars and epochs')
        while maxnan > 0:
            maxnan0 = 0
            delvec0 = None
            for i in range(epochs):
                objidvec = M[i, :]
                nancount0 = 1.0 * objidvec[np.isnan(objidvec)].size / len(objidvec)
                if nancount0 > maxnan0:
                    maxnan0 = nancount0
                    delvec0 = i

            maxnan1 = 0
            delvec1 = None
            for j in range(numstars):
                starvec = M[:, j]
                nancount1 = 1.0 * starvec[np.isnan(starvec)].size / len(starvec)
                if nancount1 > maxnan1:
                    maxnan1 = nancount1
                    delvec1 = j

            if maxnan0 >= maxnan1 and maxnan0 > 0:
                M = np.delete(M, delvec0, 0)
                self.objids.pop(delvec0)
                self.hjds = np.delete(self.hjds, delvec0)
                print('e%3d, ' % delvec0, end=' ')
            elif maxnan1 > maxnan0:
                M = np.delete(M, delvec1, 1)
                self.stars.pop(delvec1)
                self.starids.pop(delvec1)
                print('s%4d, ' % delvec1, end=' ')
            else:
                print('done (%.4f, %.4f)' % (maxnan0, maxnan1))

            maxnan = max(maxnan0, maxnan1)
            epochs, numstars = M.shape

        self._savematrix(M, 'reduced')

        print('saved reduced matrix', M.shape)

        self._saveimage(self.filename + '_reducedmatrix.png', M, sort='mean')

    def clean(self, twosigma=False, sigmaclip=5.0):
        """
        cleans the reduced matrix by performing the following tasks:
        * subtract the mean magnitude for each star --> O
        * subtract the mean magnitude for each epoch --> P
        * remove the epochs where the stddev is larger than 2 sigmas
        * calculate the PCA and remove the first component
        * subtract the mean magnitude for each epoch where the stddev for the 
          star is below 0.01 mag --> Q
        """

        M = self._loadmatrix('reduced')

        print('matrix to clean:', M.shape)

        # calculate the mean for each lightcurve-vector
        meanvec = np.mean(M, axis=0)

        O = np.ones(M.shape)
        O = O * meanvec

        # calculate the mean for each epoch
        meanvec1 = np.mean(M - O, axis=1)

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

            m1std = np.std(M1) * 2.0
            stdvec = np.where(np.std(M1, axis=1) > m1std)[0]

            fig, [ax1, ax2] = plt.subplots(1, 2)
            ax1.plot(np.sort(np.std(M1, axis=0)))
            ax1.axhline(m1std)
            ax1.set_ylabel('std')
            ax1.set_xlabel('stars')

            ax2.plot(np.sort(np.std(M1, axis=1)))
            ax2.axhline(m1std)
            ax2.set_xlabel('epochs')
            fig.savefig(self.filename + '_matrixstat.pdf')
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
                print(stdvec)

            print('matrix after 2sigma epoch reduction:', M1.shape)
        M = M1
        refstar = np.argmin(np.std(M, axis=0))
        print('reference star: %d, std = %.4f' % (refstar, np.std(M[:, refstar])))

        # perform sigma clipping at a level of 3 sigmas
        # m1std = np.std(M1)

        # i = np.where(M1 > sigmaclip*m1std)
        # M1[i] = sigmaclip*m1std
        # j = np.where(M1 < -sigmaclip*m1std)
        # M1[j] = -sigmaclip*m1std

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
        # Mhat = np.dot(U, np.dot(S, V.T))
        # print "Using all PCs, MSE = %.6G" %(np.mean((M - Mhat)**2))

        # if we use only the first n PCs the reconstruction is less accurate
        level = 1
        Mhat2 = np.dot(U[:, :level], np.dot(S[:level, :level], V[:, :level].T))
        # print "Using first %d PCs, MSE = %.6G" %(level, np.mean((M - Mhat2)**2))
        # M = M - Mhat2
        M = np.dot(U[:, level:], np.dot(S[level:, level:], V[:, level:].T))

        # calculate where the std for each lightcurve-vector is < 0.01 mag
        stdvec = np.std(M, axis=0)
        i = np.where(stdvec < 0.01)[0]
        print('stars with std<0.01: %d' % len(i))
        # calculate the mean for each epoch
        meanvec2 = np.mean(M[:, i], axis=1)

        Q = np.ones(M.shape)
        Q = (Q.T * meanvec2).T

        M = M - Q
        refstar = np.argmin(np.std(M, axis=0))
        print('new reference star: %d, std = %.4f' % (refstar, np.std(M[:, refstar])))

        # test for trends in matrix: ##########################################

        #         i = np.argsort(meanvec)
        #         R = M[:, i]
        #         for i in range(epochs):
        #             #plt.scatter(np.arange(numstars), R[i, :], edgecolor='none')
        #             x = np.arange(numstars)
        #             z = np.polyfit(x, R[i, :], 1)
        #             x = [0.0, numstars]
        #             #plt.plot(x, np.polyval(z, x))
        #             R[i, :] -= np.polyval(z, np.arange(numstars))
        #             #plt.show()
        #         plt.imshow(R)
        #         plt.show()
        #         self._saveimage(self.filename+'_cleanedmatrixR.png', R, sort='std')

        # ######################################################################

        self._savematrix(M, 'cleaned')

        plt.semilogy(s[0:30])
        plt.title('singular values')
        plt.savefig(self.filename + '_cleanedmatrixs.pdf')
        plt.close()

        self._saveimage(self.filename + '_cleanedmatrix.png', M, sort='std')
        self._saveimage(self.filename + '_cleanedmatrixM1.png', M1, sort='std')
        self._saveimage(self.filename + '_cleanedmatrixMhat.png', Mhat2, sort='std')

    def save_lightcurves(self):
        """
        save the lightcurves for all starids
        """
        M = self._loadmatrix('cleaned')
        epochs, _ = M.shape

        a = np.empty([epochs, 2])

        for starid in self.starids:
            print('.', end=' ')
            i = self.starids.index(starid)
            filename = self.lightcurvepath + '/%s.dat' % starid
            a[:, 0] = self.hjds
            a[:, 1] = M[:, i]
            np.savetxt(filename, a, fmt='%.5f %.4f')
        print('done')

    def _plot_lightcurve(self, starid, t, m, axis):
        """
        plot the lightcurve for a given star
        """
        m -= np.nanmean(m)
        t -= min(t)
        i = np.where(np.isfinite(m))
        t = t[i]
        m = m[i]
        plt.xticks(np.arange(0, 200, 10))
        plt.yticks(np.arange(-0.5, 0.5, 0.01))
        plt.xlim(min(t), max(t))
        plt.scatter(t, m, edgecolor='none', facecolor='g', s=5)
        plt.plot(t, m, 'gray')
        std = 5.0 * np.std(m)
        plt.text(0.95, 0.95, '#' + starid.split('#')[1],
                 fontsize=12,
                 verticalalignment='top',
                 horizontalalignment='right',
                 transform=axis.transAxes)
        plt.ylim([std, -std])

    def make_lightcurves(self, show=False):
        """plot lightcurve"""

        fig_width = 18.3 / 2.54  # width in inches, was 7.48in
        fig_height = 23.3 / 2.54  # height in inches, was 25.5
        fig_size = [fig_width, fig_height]
        # set plot attributes
        params = {'backend': 'Agg',
                  'axes.labelsize': 8,
                  'axes.titlesize': 10,
                  'font.size': 8,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'figure.figsize': fig_size,
                  'savefig.dpi': 300,
                  'font.family': 'sans-serif',
                  'axes.linewidth': 0.5,
                  'xtick.major.size': 2,
                  'ytick.major.size': 2,
                  }
        rcParams.update(params)

        M = self._loadmatrix('cleaned')

        sp = 1
        lc = 1
        rows = 7
        cols = 4
        for starid in self.starids:
            i = self.starids.index(starid)
            ax = plt.subplot(rows, cols, sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            if len(np.where(np.isfinite(M[:, i]))) == 0:
                continue
            self._plot_lightcurve(starid, self.hjds, M[:, i], ax)
            if sp == rows * cols + 1 or starid == self.starids[-1]:
                plt.tight_layout()
                if show:
                    plt.show()
                else:
                    print(self.plotpath + self.field + '_lc%d.pdf' % lc)
                    plt.savefig(self.plotpath + self.field + '_lc%d.pdf' % lc)
                lc += 1
                sp = 1
                plt.close()
