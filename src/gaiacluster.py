#!/usr/bin/env python
# -*- coding: utf-8 -*-
from cluster import Cluster
#from matplotlib.mlab import dist

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2018 organization_name"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__="2018-06-07"

from matplotlib import pyplot as plt
import numpy as np

class GaiaCluster(object):
    '''
    classdocs
    '''


    def __init__(self, clustername, ra=None, dec=None, diam=None, dist=None, pmra=None, pmdec=None, pmradius=5.0):
        '''
        Constructor
        '''
        self.clustername = clustername
        self.cluster  = Cluster(clustername)
        self.coordinates = self.cluster.coordinates
        
        print self.cluster
        
        def setnotNone(default, clustervalue):
            if default is None:
                if clustervalue is None:
                    raise ValueError
                return clustervalue
            else:
                return default
        
        self.ra = setnotNone(ra, self.cluster['ra'])
        self.dec = setnotNone(dec, self.cluster['dec'])
        self.diam = setnotNone(diam, self.cluster['diam'])
        self.dist = setnotNone(dist, self.cluster['d'])
        self.pmra = setnotNone(pmra, self.cluster['pmra'])
        self.pmdec = setnotNone(pmdec, self.cluster['pmdec'])
        
        self.pmradius = pmradius
        
    def fit_gaussian(self, parameter, mean=None):
        from functions import gauss_fit
        
        hist, bin_edges = np.histogram(parameter)
        bin_centers = ((np.roll(bin_edges, -1) + bin_edges) / 2.0)[:-1]
        
        return gauss_fit(bin_centers, hist, mean=mean)
        
        
    def gaiaquery(self):
        '''
        The is the code sample to access https://gaia.aip.de/tap with your token account. 
        Please copy this file and replace 'your_token' with your authorization token from https://gaia.aip.de/accounts/token/
        '''
        
        from daiquiri_tap import DaiquiriTap
        TOKEN = '36567049f50f830ab96cd32f89726295ddca964f'
        service = DaiquiriTap(url='https://gaia.aip.de/tap', token=TOKEN)

        #Here you create the job with the query.
        # WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 282.825, 10.3183, 0.5)) = 1
        #, phot_g_mean_mag "gmag",  phot_bp_mean_mag "bmag", phot_rp_mean_mag "rmag"
        #     AND (brmag > -0.5) AND (brmag < 2.0)
        params = {'name': self.clustername,
                  'ra': self.ra,
                  'dec': self.dec,
                  'diam': self.diam/60.0,
                  'radius': self.diam/120.0,
                  'dist': self.dist,
                  'pmra': self.pmra,
                  'pmdec': self.pmdec,
                  'pmradius': self.pmradius}
        querystring = '''SELECT TOP 10000 source_id, ra, dec, phot_g_mean_mag AS gmag, phot_bp_mean_mag-phot_rp_mean_mag AS brmag, pmra, pmdec, 1000./parallax AS dist
        FROM gdr2.gaia_source 
        WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %(ra)f, %(dec)f, %(radius)f)) = 1 
        AND phot_g_mean_mag < 18.5 
        AND phot_bp_mean_mag-phot_rp_mean_mag > 0.0 AND phot_bp_mean_mag-phot_rp_mean_mag < 2.5
        AND 1000./parallax > %(dist)f*0.85 AND 1000./parallax < %(dist)f*1.15
        AND pmra > %(pmra)f-%(pmradius)f AND pmra <= %(pmra)f+%(pmradius)f
        AND pmdec > %(pmdec)f-%(pmradius)f AND pmdec <= %(pmdec)f+%(pmradius)f''' % params
        #print querystring
        job = service.launch_job(querystring)
        #print(job)
        r = job.get_results()
        #newdist = np.mean(r['dist'])  
        
        _, newra, newra_sigma = self.fit_gaussian(r['ra'], mean=self.ra)
        _, newdec, newdec_sigma = self.fit_gaussian(r['dec'], mean=self.dec)
        _, newdist, newdist_sigma = self.fit_gaussian(r['dist'], mean=self.dist)
        _, newpmra, newpmra_sigma = self.fit_gaussian(r['pmra'], mean=self.pmra)
        _, newpmdec, newpmdec_sigma = self.fit_gaussian(r['pmdec'], mean=self.pmdec)
                                
 
        plt.style.use('a4paper.mplstyle')
        #newpmra = np.median(r['pmra']) 
        #newpmdec = np.median(r['pmdec']) 
        plt.subplot(221)
        plt.title(self.clustername)
        plt.plot(r['brmag'], r['gmag'], 'k.')
        plt.xlabel('b - r')
        plt.ylabel('G mag')
        plt.ylim(plt.ylim()[::-1])
        
        ax2 = plt.subplot(222)
        ax2.set_aspect(1.)
        
        plt.plot(r['ra'], r['dec'],'k.')
        plt.xlabel('RA %.2f +- %.2f' % (newra, newra_sigma))
        plt.ylabel('Dec %.2f +- %.2f'% (newdec,newdec_sigma))
        plt.axvline(self.cluster['ra'], color='r', linestyle='--')
        plt.axhline(self.cluster['dec'], color='r', linestyle='--')
        plt.axvline(newra, color='g', linestyle='--')
        plt.axhline(newdec, color='g', linestyle='--')
        
        plt.xlim(plt.xlim()[::-1])
        
        ax3 = plt.subplot(223)
        ax3.set_aspect(1.)
        plt.plot(r['pmra'],r['pmdec'], 'k.', alpha=0.5)
        plt.axvline(self.cluster['pmra'], color='r', linestyle='--')
        plt.axhline(self.cluster['pmdec'], color='r', linestyle='--')
        plt.axvline(newpmra, color='g', linestyle='--')
        plt.axhline(newpmdec, color='g', linestyle='--')
        plt.xlim(self.pmra-self.pmradius, self.pmra+self.pmradius)
        plt.ylim(self.pmdec-self.pmradius, self.pmdec+self.pmradius)
        plt.xlabel('pm RA %.2f +- %.2f' % (newpmra, newpmra_sigma))
        plt.ylabel('pm Dec %.2f +- %.2f'% (newpmdec,newpmdec_sigma))
        plt.grid()
        
        plt.subplot(224)
        plt.title('total: %d' % len(r['source_id']))
        plt.hist(r['dist'], bins=20)
        plt.axvline(self.cluster['d'], color='r', linestyle='-.')
        plt.axvline(newdist, color='g', linestyle='-.')
        
        plt.xlabel('distance %.1f +- %.1f' % (newdist, newdist_sigma))
        plt.savefig('/work2/jwe/SOCS/plots/Gaia DR2 %s.pdf' % self.clustername)
        plt.close()
        print("GaiaCluster('%s', ra=%.2f, dec=%.2f, pmra=%.2f, pmdec=%.2f, diam=120.0, pmradius=2.0, dist=%.1f)" % \
              (self.clustername, newra, newdec, newpmra, newpmdec, newdist))
        
        
if __name__ == '__main__':
    #c = GaiaCluster('NGC 6633', ra=276.86, dec=6.66, dist=395.3, pmra=1.18, pmdec=-1.74, diam=180.0, pmradius=1.0)
    #c.gaiaquery()
    #c = GaiaCluster('NGC 6709', ra=282.84, dec=10.34, pmra=1.45, pmdec=-3.52, dist=1090.1, pmradius=1.0, diam=90.0)
    #c.gaiaquery()
    #c = GaiaCluster('NGC 2323', ra=105.69, dec=-8.36, pmra=-0.76, pmdec=-0.64, dist=1003.7, pmradius=2.0, diam=45.0)
    #c.gaiaquery()
    #c = GaiaCluster('NGC 6940', diam=120.0, pmradius=8.0)
    #c = GaiaCluster('IC 4756', pmra=1.26, pmdec=-4.90, pmradius=2.0, diam=120.0, dist=478.5)
    #c = GaiaCluster('NGC 2422', ra=114.14, dec=-14.47, pmra=-7.05, pmdec=1.00, diam=120.0, pmradius=2.0, dist=482.2)
    #c = GaiaCluster('NGC 2236', ra=97.41, dec=6.84, pmra=-0.77, pmdec=0.05, diam=20.0, pmradius=1.1, dist=2749.2)
    #c = GaiaCluster('NGC 1647', ra=71.50, dec=19.12, pmra=-1.01, pmdec=-1.53, diam=120.0, pmradius=2.0, dist=593.9)
    #c = GaiaCluster('NGC 2281', ra=102.09, dec=41.04, pmra=-2.95, pmdec=-8.31, diam=120.0, pmradius=2.0, dist=529.0)
    c = GaiaCluster('NGC 1528', ra=63.87, dec=51.21, pmra=2.22, pmdec=-2.25, diam=90.0, pmradius=1.5, dist=1054.2)
    c.gaiaquery()
    
    