#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 9, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
import numpy as np
from tools import log
class Calibration(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.data = []
        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost, dictcursor=True)
 
    def loaddb(self):
        query = """SELECT * 
        FROM referencestars 
        WHERE starid LIKE 'M67%%' 
        ORDER BY starid"""
        self.data = self.wifsip.query(query)

    def getphot(self, ra, dec, filtercol):
        from datasource import DataSource
        wifsip = DataSource(database=config.dbname, 
                                 user=config.dbuser, 
                                 host=config.dbhost,
                                 dictcursor=True)
        params = {'filtercol': filtercol,
                  'ra': ra,
                  'dec': dec}
        #'IC 4756 v2 %% uvby' 
        query = """SELECT mag_auto, zeropnt, expt, flux_auto
        FROM phot, frames
        WHERE object like 'M 67%%'
        AND phot.objid = frames.objid
        AND filter='%(filtercol)s'
        AND circle(phot.coord,0.0)<@circle(point(%(ra).11f,%(dec).11f),0.6/3600.0)
        AND flags=0""" % params
        data = wifsip.query(query) 
        zeropnt = np.array([d['zeropnt'] for d in data])
        mag_isocor = np.array([d['mag_auto'] for d in data])
        expt = np.array([d['expt'] for d in data])
        flux_auto = np.array([d['flux_auto'] for d in data])
        if len(mag_isocor)==0:
            return np.nan, np.nan, np.nan
        #mags = mag_isocor-zeropnt
        mags = -2.5*np.log10(flux_auto/expt)
        valid = abs(mags - np.mean(mags)) < 2.0*np.std(mags)
        mags = np.compress(valid, mags)
        mag = np.mean(mags)

        return mag, np.std(mags), len(mags)

        
    def calstar(self, star):
        #print '%s %.6f %.6f' % (star['starid'],star['ra'], star['dec']),
        mag = {}
        err = {}
        n = {}
        for filtercolor in ['b','y']:
            mag[filtercolor], err[filtercolor], n[filtercolor] = self.getphot(star['ra'], star['dec'], filtercolor)
        if n['b']>0 and n['y']>0:
            try:
                log(config.datapath+'cal.dat', '%s %.6f %.6f %.3f %.4f %.4f %3d %.4f %.4f %3d' % \
                    (star['starid'],star['ra'], star['dec'], star['b-y'], mag['b'],err['b'],n['b'],mag['y'],err['y'],n['y']))
            except TypeError:
                print star['starid']
            
    
    def calstars(self):
        log(config.datapath+'cal.dat', '#starid RA Dec by bmag berr nb ymag yerr ny')
        for star in self.data:
            self.calstar(star)
 
if __name__ == '__main__':
    cal = Calibration()
    cal.loaddb()
    cal.calstars()