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
'''
select name, ra, dec, diam, d, ebv, logage, pmra, pmdec, nc from clusters where dec>-20 and ebv<0.3 and d<2000 and logage>8.0 and logage< 9.0 and nc>100 order by name limit 100;
select name,ra,dec,ebv, d, logage, nc from clusters where ebv<0.31 and d<1600 and logage>8.0 and logage< 9.5 and nc>25 order by name;
   name   |  ebv  |  d   | logage |  nc  
----------+-------+------+--------+------
 NGC 752  | 0.034 |  457 |   9.05 |   50
 NGC 1528 |  0.26 | 1090 |    8.6 |  405
 NGC 1647 |  0.37 |  540 |  8.158 |  656
 NGC 2126 |  0.27 | 1090 |    9.1 |   68
 NGC 2236 | 0.479 | 2930 |  8.538 |  184
 NGC 2301 |  0.03 |  870 |    8.2 |  355
 NGC 2281 | 0.063 |  558 |  8.554 |  330
 NGC 2324 |  0.25 | 3800 |   8.65 |  316
 NGC 6709 | 0.304 | 1075 |  8.178 |  700
 NGC 6866 |   0.1 | 1470 |    8.8 |  498
 IC 4756  | 0.192 |  484 |  8.699 |   30
 NGC 6633 | 0.182 |  376 |  8.629 |  872
 NGC 1912 |  0.25 | 1400 |    8.5 | 2579
 NGC 2437 |   0.1 | 1510 |    8.4 | 1287
 NGC 6940 | 0.214 |  770 |  8.858 |  877
 NGC 2422 |  0.07 |  490 |  7.861 |   23
 NGC 2548 |  0.03 |  770 |    8.6 |  522
 NGC 2323 |   0.2 |  950 |      8 |  377
 NGC 2682 |  0.03 |  808 |   9.45 |  399

'''
class GaiaCluster(object):
    '''
    Gaia Cluster class
    '''


    def __init__(self, clustername, ra=None, dec=None, diam=None, dist=None, pmra=None, pmdec=None, pmradius=5.0):
        '''
        Constructor
        '''
        self.clustername = clustername
        self.cluster  = Cluster(clustername)
        self.coordinates = self.cluster.coordinates
        
        print(self.cluster)
        
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
        from _config import TOKEN
        service = DaiquiriTap(url='https://gaia.aip.de/tap', token=TOKEN)

        #Here you create the job with the query.
        # WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 282.825, 10.3183, 0.5)) = 1
        #AND pmra > %(pmra)f-%(pmradius)f AND pmra <= %(pmra)f+%(pmradius)f
        #AND pmdec > %(pmdec)f-%(pmradius)f AND pmdec <= %(pmdec)f+%(pmradius)f
        params = {'name': self.clustername,
                  'ra': self.ra,
                  'dec': self.dec,
                  'diam': self.diam/60.0,
                  'radius': self.diam/120.0,
                  'dist': self.dist,
                  'pmra': self.pmra,
                  'pmdec': self.pmdec,
                  'pmradius': self.pmradius}
        querystring = '''SELECT TOP 100000 source_id, ra, dec, phot_g_mean_mag AS gmag, bp_rp AS brmag, pmra, pmdec, 1000./parallax AS dist
        FROM gdr2.gaia_source 
        WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %(ra)f, %(dec)f, %(radius)f)) = 1 
        AND phot_g_mean_mag < 18.5 
        AND phot_bp_mean_mag-phot_rp_mean_mag > 0.0 AND phot_bp_mean_mag-phot_rp_mean_mag < 2.5
        AND 1000./parallax > %(dist)f*0.5 AND 1000./parallax < %(dist)f*1.5
        AND CONTAINS(POINT('', pmra, pmdec), CIRCLE('', %(pmra)f, %(pmdec)f, %(pmradius)f)) = 1''' % params
        #name = 'Gaia DR2 SOCS %s' % self.clustername
        output_file = '/work2/jwe/SOCS/data/%s.votable' % self.clustername
        job = service.launch_job(querystring, output_file=output_file, dump_to_file=True)
        #print(job)
        self.data = job.get_results()
        #job.set_output_file(output_file)
        #job.save_results(verbose=True)
        
        
    
    def plot(self):
        r = self.data
        _, newra, newra_sigma = self.fit_gaussian(r['ra'], mean=self.ra)
        _, newdec, newdec_sigma = self.fit_gaussian(r['dec'], mean=self.dec)
        _, newdist, newdist_sigma = self.fit_gaussian(r['dist'], mean=self.dist)
        _, newpmra, newpmra_sigma = self.fit_gaussian(r['pmra'], mean=self.pmra)
        _, newpmdec, newpmdec_sigma = self.fit_gaussian(r['pmdec'], mean=self.pmdec)
        newdiam = np.max([newra_sigma, newdec_sigma]) * 3.0 * 60.0  
        newpmradius =  np.max([newpmra_sigma, newpmdec_sigma]) * 3.0 * 2.0          

        plt.style.use('a4paper.mplstyle')
        #newpmra = np.median(r['pmra']) 
        #newpmdec = np.median(r['pmdec']) 
        plt.subplot(221)
        plt.title(self.clustername)
        plt.plot(r['brmag'], r['gmag'], 'k.', alpha=0.25, mec='None')
        plt.xlabel('b - r')
        plt.ylabel('G mag')
        plt.ylim(plt.ylim()[::-1])
        plt.minorticks_on()
        plt.grid()
        
        ax2 = plt.subplot(222)
        ax2.set_aspect(1.)
        
        plt.plot(r['ra'], r['dec'],'k.', alpha=0.5, mec='None')
        plt.xlabel('RA %.2f +- %.2f' % (newra, newra_sigma))
        plt.ylabel('Dec %.2f +- %.2f'% (newdec,newdec_sigma))
        plt.axvline(self.cluster['ra'], color='r', linestyle='--')
        plt.axhline(self.cluster['dec'], color='r', linestyle='--')
        plt.axvline(newra, color='g', linestyle='--')
        plt.axhline(newdec, color='g', linestyle='--')
        plt.grid()
        plt.xlim(plt.xlim()[::-1])
        
        ax3 = plt.subplot(223)
        ax3.set_aspect(1.)
        plt.plot(r['pmra'],r['pmdec'], 'k.', alpha=0.5, mec='None')
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
        plt.savefig('/work2/jwe/SOCS/plots/Gaia DR2 new %s.pdf' % self.clustername)
        plt.close()
        print("GaiaCluster('%s', ra=%.2f, dec=%.2f, pmra=%.2f, pmdec=%.2f, diam=%.1f, pmradius=%.1f, dist=%.1f)" % \
              (self.clustername, newra, newdec, newpmra, newpmdec, newdiam, newpmradius, newdist))
        
        
if __name__ == '__main__':
    #c = GaiaCluster('NGC 752', ra=29.28, dec=37.77, pmra=9.86, pmdec=-11.69, diam=180.0, pmradius=2.0, dist=448.9)
    #c = GaiaCluster('NGC 1528', ra=63.87, dec=51.21, pmra=2.22, pmdec=-2.25, diam=90.0, pmradius=1.1, dist=1054.2)
    #c = GaiaCluster('NGC 1647', ra=71.50, dec=19.12, pmra=-1.01, pmdec=-1.53, diam=120.0, pmradius=2.0, dist=593.9)
    #c = GaiaCluster('NGC 1912', pmra=1.49, pmdec=-4.48, diam=120.0, pmradius=2.0) # double cluster
    c = GaiaCluster('NGC 1996', pmra=0.94, pmdec=-3.21, dist=1710.5 )
    #c = GaiaCluster('NGC 2126', pmra=1.49, pmdec=-3, diam=120.0, pmradius=4.0)
    #c = GaiaCluster('NGC 2236', ra=97.41, dec=6.84, pmra=-0.77, pmdec=0.05, diam=20.0, pmradius=1.1, dist=2749.2)
    #c = GaiaCluster('NGC 2281', ra=102.09, dec=41.04, pmra=-2.95, pmdec=-8.31, diam=120.0, pmradius=2.0, dist=529.0)
    #c = GaiaCluster('NGC 2301', ra=102.94, dec=0.46, pmra=-1.35, pmdec=-2.16, diam=120.0, pmradius=2.0, dist=877.3)
    #c = GaiaCluster('NGC 2323', ra=105.69, dec=-8.36, pmra=-0.76, pmdec=-0.64, dist=1003.7, pmradius=2.0, diam=45.0)
    #c = GaiaCluster('NGC 2324', ra=106.03, dec=1.03, pmra=-0.27, pmdec=-0.09, diam=30.0, pmradius=2.0, dist=3716.1)
    #c = GaiaCluster('NGC 2422', ra=114.14, dec=-14.47, pmra=-7.05, pmdec=1.00, diam=120.0, pmradius=2.0, dist=482.2)
    #c = GaiaCluster('NGC 2437', ra=115.46, dec=-14.84, pmra=-3.81, pmdec=0.41, diam=90.0, pmradius=1.0, dist=1633.1) # M 46
    #c = GaiaCluster('NGC 2548', ra=123.40, dec=-5.73, pmra=-1.31, pmdec=1.03, diam=120.0, pmradius=2.0, dist=777.6) # M 48
    #c = GaiaCluster('NGC 2682', ra=132.85, dec=11.82, pmra=-10.98, pmdec=-2.99, diam=120.0, pmradius=2.0, dist=879.4) # M 67
    #c = GaiaCluster('IC 4756', ra=279.69, dec=5.42, pmra=1.26, pmdec=-4.90, diam=74.8, pmradius=1.9, dist=478.3)
    #c = GaiaCluster('NGC 6633', ra=276.86, dec=6.66, dist=395.3, pmra=1.18, pmdec=-1.74, diam=180.0, pmradius=2.0)
    #c = GaiaCluster('NGC 6709', ra=282.84, dec=10.34, pmra=1.45, pmdec=-3.52, dist=1090.1, pmradius=1.0, diam=90.0)
    #c = GaiaCluster('NGC 6866', ra=300.99, dec=44.16, pmra=-1.36, pmdec=-5.81, diam=70, pmradius=1.0, dist=1455.5)
    #c = GaiaCluster('NGC 6940', ra=308.64, dec=28.27, pmra=-1.97, pmdec=-9.42, diam=63.3, pmradius=2.0, dist=1051.3) #20 34 26.0 +28 17 00
    
    #c = GaiaCluster('NGC 2168', ra=92.26, dec=24.34, pmra=2.33, pmdec=-2.88, diam=32.3, pmradius=2.8, dist=882.5)
    #c = GaiaCluster('NGC 1039', ra=40.54, dec=42.74, pmra=0.70, pmdec=-5.72, diam=34.8, pmradius=2.4, dist=510.1)
    #c = GaiaCluster('NGC 1996', ra=84.57, dec=25.81, pmra=1.08, pmdec=-3.29, diam=39.1, pmradius=17.6, dist=1487.1)
    
    #c = GaiaCluster('NGC 2353', ra=108.63, dec=-10.26, pmra=-0.94, pmdec=0.82, diam=37.4, pmradius=2.0, dist=1227.3)
    #c = GaiaCluster('NGC 2374', ra=110.98, dec=-13.23, pmra=-4.6, pmdec=0.62, diam=72.3, pmradius=1.1, dist=1333.0)
    #c = GaiaCluster('NGC 2539', ra=122.66, dec=-12.83, pmra=-2.50, pmdec=-0.65, diam=40.6, pmradius=2.0, dist=1327.8)
    #c = GaiaCluster('NGC 6811', ra=294.35, dec=46.38, pmra=-3.22, pmdec=-8.92, diam=57.5, pmradius=2.0, dist=1151.9)
    #c = GaiaCluster('NGC 7082', ra=322.25, dec=47.12, pmra=-0.41, pmdec=-0.99, diam=36.0, pmradius=3.2, dist=1426.9)
    #c = GaiaCluster('NGC 7209', ra=331.24, dec=46.49, pmra=2.41, pmdec=0.21, diam=77.8, pmradius=2.0, dist=1216.9)
    #c = GaiaCluster('NGC 7243', ra=333.75, dec=49.84, pmra=0.46, pmdec=-2.81, diam=71.2, pmradius=1.7, dist=888.3)
    
    #c = GaiaCluster('NGC 2126', ra=90.62, dec=49.89, pmra=0.96, pmdec=-2.66, diam=30, pmradius=1.5, dist=1318.7)
    #c = GaiaCluster('NGC 2184', diam=120.0, dist=995.6)
    c.gaiaquery()
    c.plot()
    #for name in ['NGC 2180','NGC 2184','NGC 2251','NGC 2252','NGC 2270','NGC 2306','NGC 2319','NGC 2353','NGC 2358','NGC 2374','NGC 2413','NGC 2423','NGC 2430','NGC 2432','NGC 2539','NGC 6596','NGC 6639','NGC 6728','NGC 6793','NGC 6811','NGC 6828','NGC 7058','NGC 7082','NGC 7209','NGC 7243']:
        #c = GaiaCluster(name, diam=60.0)
        #c.gaiaquery()
    
    