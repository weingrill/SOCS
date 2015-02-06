#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 21, 2013

@author: Joerg Weingrill <jweingrill@aip.de>

Data reduction Class for NGC2236
'''

import config
import logging
from ngc2236star import NGC2236Star            
import pylab as plt
import numpy as np
logging.basicConfig(filename=config.projectpath+'m48_analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')


class Ngc2236(object):
    def __init__(self):
        """Constructor"""
        from cluster import Cluster
        from astronomy import mag_distance
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.stars = []
         
        c = Cluster('NGC 2236')
        self.age = 10**c['logage']/1e6 # in Myr
        self.ebv = c['ebv']
        self.dm = mag_distance(c['d'])- 4.83

    def clearperiods(self):
        """
        reset the periods in the database table
        """
        if not raw_input('press Y to erase the periods in table ')=='Y':
            return
        query="""UPDATE ngc2236 
        SET period=NULL, period_err=NULL, amp=NULL, amp_err=NULL
        WHERE period>0;
        """
        logger.info('resetting periods ...')
        self.wifsip.execute(query)

    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        query = """SELECT starid 
        FROM ngc2236 
        WHERE NOT corotid IS NULL;"""
        
        logger.info('fetching stars ...')
        result = self.wifsip.query(query)
        logger.info('... %d stars found' % len(result))
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]
    
    def analysis(self, show=False):
        """perform a PDM analysis on each lightcurve"""
        from matplotlib import rcParams
        from functions import sigma_clip, phase
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

        for starid in self.stars:
            star = NGC2236Star(starid)
            print '%-24s '% starid,
            lc = star.lightcurve()
            meanflux = np.nanmean(lc.flux)
            lc.flux /= meanflux
            lc.flux *= meanflux
            lc.rebin(0.0125)
            lc.interpolate()
            lc.normalize()
            raw_time, raw_flux = (lc.time, lc.flux)
            
            lc.medfilter()
            lc.jmpflt()
            lc.detrend()
            time, flux = (lc.time, lc.flux)
            time -= min(time)
            
            # convert to magnitudes
            mag = -2.512*np.log10(meanflux*flux) + 24
            # remove the mean so we are around 0.0 magnitudes
            mag -= np.mean(mag)
            
            # perform a 3sigma clipping
            time, mag = sigma_clip(time, mag)

            # calculate fourier spectrum with zero padding
            n = len(mag)
            n2 = 8*n
            ft = np.fft.rfft(mag, n2)
            #amp = abs(ft)/n
            power = abs(ft)**2/n
            
            freq = np.fft.fftfreq(n2, d=lc.dt())
            
            
            # apply a bessel filter to eliminate runlength artifacts
            from scipy import signal
            
            b, a = signal.butter(4, 2./max(time), 'high', analog=True)
            _, h = signal.freqs(b, a, worN=freq[1:n])
            filt_pwr = abs(ft[1:n]*h)**2/n

            i = np.argmax(filt_pwr)+1
            maxfreq = freq[i]
            period = 1./maxfreq
            
            from pdm import pdm

            #refine period using pdm
            pdm_periods, pdm_thetas = pdm(time, mag, period/1.5, period*1.5, 0.0125)
            
            i = np.argmin(pdm_thetas)
            period = pdm_periods[i]
            
            print 'P = %5.2f' % period,
            star['period'] = period
            periods = 1./freq[1:n]
            periods = periods[::-1]
            filt_pwr = filt_pwr[::-1]
            
            norm = np.mean(filt_pwr)
            print '%.1f' % (max(filt_pwr)/norm)
            power[1:n] /= norm
            filt_pwr /= norm
            star['amp'] = max(filt_pwr)/norm
            ph_time, ph_mag = phase(time, mag, period)
            num = len(ph_time)/round(max(time)/period)
            rph_mag, _ = signal.resample(ph_mag, num, t = ph_time)
            rph_time =  np.linspace(min(ph_time),max(ph_time),num)
            bv = star['bv']
            if bv is None: bv=-99
            
            plt.subplot(411) ##################################################
            plt.title('%s (%d) B-V=%.2f' % (starid, star['corotid'], bv))
            plt.scatter(raw_time-min(raw_time), raw_flux, edgecolor='none', alpha=0.5, s=3, color='k')
            plt.ylim(min(raw_flux),max(raw_flux))
            plt.xlim(0,max(time))
            
            plt.subplot(412) ##################################################
            plt.plot(time, -mag, 'k')
            plt.ylim(min(-mag),max(-mag))
            plt.xlim(0,max(time))
            
            plt.subplot(413) ##################################################
            plt.plot(1./freq[1:n], power[1:n], 'k')
            plt.plot(1./freq[1:n], filt_pwr, 'g')
            #plt.plot(1./w, 20 * np.log10(abs(h)))
            plt.axvline(period)
            #plt.axhline(np.mean(filt_pwr[1:n]))
            plt.xlim(0.1, max(time))
            
            plt.subplot(414) ##################################################
            plt.plot(rph_time, -rph_mag,color='k')
            plt.plot(rph_time+period, -rph_mag,color='k')
            #plt.plot(phased_time, phased_mag, 'g')
            #plt.plot(phased_time+period, phased_mag, 'g')
            plt.axvline(period, linestyle='--')
            plt.xlabel('P = %.2f' % period)
            plt.xlim(0., period*2)
            
            #plt.grid()

            if show:
                plt.show()
            else:
                plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['corotid']))
            plt.close()

    def set_tab_column(self):
        
        self.wifsip.execute('UPDATE ngc2236 set tab=NULL;')
        query = """SELECT starid
            FROM ngc2236 
            WHERE NOT bv IS NULL
            ORDER BY vmag;"""
        result = self.wifsip.query(query)
        starids = [r[0] for r in result]
        
        for starid in starids:
            tab = starids.index(starid)+1
            print '%4d %s' % (tab,starid)
            query = "UPDATE ngc2236 set tab=%d WHERE starid='%s';" % (tab,starid)
            self.wifsip.execute(query)
    
    def update_coordinates(self):
        from datasource import DataSource
        
        corot = DataSource(database='corot', user='sro', host='pina.aip.de')
        
        query = """SELECT corotid
            FROM ngc2236 
            WHERE NOT corotid IS NULL
            ORDER BY vmag;"""
        corotids = [c[0] for c in self.wifsip.query(query)]
        
        for corotid in corotids:
            query = """SELECT alpha, delta
            FROM corot 
            WHERE corotid = %d;""" % corotid
        
            ra, dec = corot.query(query)[0]
            print corotid,ra,dec
            record = dict(zip(['corotid','ra','dec'], [corotid,ra,dec]))
            query = """UPDATE ngc2236 
            SET ra = %(ra)f, dec = %(dec)f, coord = point(%(ra)f,%(dec)f)
            WHERE corotid = %(corotid)d;""" % record
            self.wifsip.execute(query, commit=False)
        self.wifsip.commit()
            
            
def calibrate(objname, filtercol):
    from calibrate import Calibrate

    calv = Calibrate(obj=objname, filtercol=filtercol)
    calv.resetframes()
    calv.getrefframes()
    for r in calv.ref: print r
    
    for ref in calv.ref:
        print 'reference:', ref
        calv.getframes(ref)
        calv.corrframes(ref)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NGC 2236 analysis')
    parser.add_argument('--calibrate', action='store_true', help='calibrate the photometry')
    parser.add_argument('--create', action='store_true', help='create table')
    parser.add_argument('--clear', action='store_true', help='clear table')
    parser.add_argument('--getframes', action='store_true', help='get frames')
    parser.add_argument('--sigmas', action='store_true', help='plot cmd')
    parser.add_argument('--bv', action='store_true', help='update B-V')
    parser.add_argument('filter', default='V', type=str, help='filter color to process')
    parser.add_argument('-a', '--analysis', action='store_true', help='analysis')
    parser.add_argument('-tab', action='store_true', help='update tab column')
    parser.add_argument('-coords', action='store_true', help='update coordinates')

    args = parser.parse_args()

    if args.calibrate: calibrate('NGC 2236 BVI', args.filter)
        
    if args.create or args.clear or args.getframes or args.sigmas or args.bv:
        from photometry import Photometry
        phot = Photometry(objname='NGC 2236 BVI', filtercol=args.filter, dbname='ngc2236')
    if args.create: phot.createtable()
    if args.clear: phot.cleartable()
    if args.getframes: phot.getframes(fields=[''])
    if args.sigmas: phot.update_sigmas()
    if args.bv: phot.update_bv()
    if args.analysis: 
        ngc2236 = Ngc2236()
        ngc2236.getstars()
        ngc2236.analysis()
        
    if args.tab: 
        ngc2236 = Ngc2236()
        ngc2236.set_tab_column()

    if args.coords: 
        ngc2236 = Ngc2236()
        ngc2236.update_coordinates()
