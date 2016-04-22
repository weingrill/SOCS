#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 19, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

from ngc6633 import LightCurve
import config
import logging
import numpy as np
#import matplotlib.pyplot as plt

logging.basicConfig(filename=config.projectpath+'rmsanalysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC6633 analysis')

class RMSAnalysis(object):
    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        self.stars = []
#        self.age = 10**8.557/1e6 # in Myr from Webda
#        self.ebv = 0.031 # from Webda
#        self.dm = 9.53 # from Webda

    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        query = """SELECT starid, vmag, bv 
        FROM ngc6633 
        WHERE (good IS NULL OR good)
        ORDER BY vmag;""" 

        result = self.wifsip.query(query)
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]
        self.vmag = [s[1] for s in result]
        self.bv = [s[2] for s in result]
    
    def setrms(self, starid, rms):
        """
        write the rms to the database
        """
        params = {'starid': starid, 'rms': rms}
        if np.isfinite(rms):
            query = """UPDATE ngc6633 SET rms = %(rms)f WHERE starid='%(starid)s';""" % params
        else: 
            query = """UPDATE ngc6633 SET rms = NULL WHERE starid='%(starid)s';""" % params
        self.wifsip.execute(query) 

    def analysis(self, show=False):
        """perform a PDM analysis on each lightcurve"""
        from matplotlib import rcParams
        print 'Analysis'

        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'font.size': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)

        
        
        for starid,vmag,bv in zip(self.stars,self.vmag,self.bv):
            print '%-24s '% starid,
            try:
                lc = LightCurve(starid)
                
            except IOError:
                logger.error("Can't load lightcurve %s" % starid)
                print 'no lightcurve'
                #self.setbad(starid)
                continue
            
            if len(lc)<50:
                logger.warn("%s: not enough datapoints" % starid)
                print 'not enough datapoints'
                continue                    
            # perform a 3sigma clipping
            lc.normalize()
            lc.detrend()
            lc.normalize()
            
            self.hjd = lc.hjd
            self.mag = lc.mag
            self.err = lc.err
            rms = np.std(self.mag)
            self.setrms(starid, rms)
            print '%.4f %.3f %.4f' % (vmag,bv,rms)


if __name__ == '__main__':
    rmsa = RMSAnalysis()
    rmsa.getstars()
    rmsa.analysis()
    