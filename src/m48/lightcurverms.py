#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 29, 2014

@author: jwe
'''

import config
import logging
from m48star import M48Star            
import pylab as plt
import numpy as np

class LightcurveRMS(object):
    '''
    Class to implement all methods for analyzing and plotting
    '''
    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')

    def load_lightcurve(self, tab=1):
        query="""SELECT mag_auto+corr
                    FROM phot, frames, m48stars
                    WHERE circle(phot.coord,0) <@ circle(m48stars.coord,1.0/3600)
                    AND tab=%d
                    AND filter='V'
                    AND frames.objid=phot.objid
                    AND phot.flags<8
                    AND NOT corr IS NULL
                    ORDER BY phot.objid;""" % tab
        result = self.wifsip.query(query)
        mags = np.array([r[0] for r in result])
        if len(mags)>=3:
            print '%4d' % tab,
            valid = abs(mags-np.mean(mags))<3.0*np.std(mags)
            magsc = np.compress(valid, mags)
            print 'mean1 = %6.3f'% np.mean(magsc),
            print 'std1 = %.3f' % np.std(magsc),
            print ' %.3f' % (np.std(magsc)/np.sqrt(len(magsc))),
            print 'n = %3d' % len(magsc) 
                 
           
            query = """INSERT INTO m48_lightcurves (tab, vmean, vstd, nv)
            VALUES (%d, %f, %f, %d);""" % (tab, np.mean(magsc), np.std(magsc), len(magsc))
            #query = """UPDATE m48_lightcurves SET bmean = %f, bstd = %f, nb = %d
            #WHERE tab=%d;""" % (np.mean(magsc), np.std(magsc), len(magsc), tab)
            query = query.replace('nan','NULL')
            self.wifsip.execute(query)
        else:
            print tab,'no data'
        if tab % 50 == 0:
            self.plot_cmd()
            
    def plot_cmd(self):
        query = """SELECT bv, vmean, vstd, good, member
            FROM m48stars, m48_lightcurves
            WHERE m48stars.tab=m48_lightcurves.tab;"""
        result = self.wifsip.query(query)
        bv = np.array([r[0] for r in result])
        vmag = np.array([r[1] for r in result])
        vstd = np.array([r[2] for r in result])
        good = [r[3] for r in result]
        member = [r[4] for r in result]
        plt.xlim(0.0,1.8)
        plt.ylim(18.0,8.0)
        plt.xlabel('B - V')
        plt.ylabel('V mag')
        plt.grid()
        s = ((vstd-np.min(vstd))/(np.max(vstd)-np.min(vstd)))*200.0
        plt.scatter(bv, vmag, s=s, color='r',alpha=0.75, label='star')
        members = np.where(member)
        goods = np.where(good)
        #print member,members,bv[members]
        plt.scatter(bv[members], vmag[members], s=s[members], color='b',alpha=0.75, label='member')
        plt.scatter(bv[goods], vmag[goods], s=s[goods], color='g',alpha=0.75, label='good')
        plt.savefig(config.plotpath+'M48_activity_cmd.pdf')
        #plt.legend()
        plt.close()
    
    @property
    def maxtab(self):
        result = self.wifsip.query('SELECT max(tab) from m48_lightcurves;')
        return result[0][0]+1

if __name__ == '__main__':
    lcrms = LightcurveRMS()
    #lcrms.plot_cmd()
    #exit()
    for tab in np.arange(lcrms.maxtab, 5121):
        lcrms.load_lightcurve(tab)