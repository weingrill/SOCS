#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 9, 2014

@author: jwe
'''
import config
import logging
import numpy as np
import matplotlib.pyplot as plt
logging.basicConfig(filename=config.projectpath+'analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC6633 analysis')

def calibrate(objname, filtercol):
    from calibrate import Calibrate

    calv = Calibrate(obj=objname+' %%', filtercol=filtercol)
    calv.resetframes()
    calv.getrefframes()
    for r in calv.ref: print r
    
    for ref in calv.ref:
        print 'reference:', ref
        calv.getframes(ref)
        calv.corrframes(ref)

class LightCurve(object):
    """
    Lightcurve object for NGC6633
    fetch lightcurve either from db or from file.
    ability to store lightcurve to file
    """
    def __init__(self, starid):
        self.starid = starid
        self.fromfile()

    def fromfile(self, filename=None):
        """
        loads the lightcurve from a file.
        if the filename is not given it is assembled from the starid
        """
        
        if filename is None:
            filename = config.lightcurvespath+self.starid+'.dat'
        logger.info('load file %s' % filename)
        self.hjd, self.mag, self.err = np.loadtxt(filename, unpack = True)
        
        logger.info('%d datapoints' % len(self.hjd))
        
        return (self.hjd, self.mag, self.err)
    
    def normalize(self):
        try:
            self.hjd -= np.min(self.hjd)
            self.mag -= np.mean(self.mag)
        except ValueError:
            return
    
    def detrend(self):
        try:
            par = np.polyfit(self.hjd, self.mag, 1)
            self.mag -= np.polyval(par, self.hjd)
        except TypeError:
            return

    def sigma_clip(self):
        from functions import sigma_clip
        self.hjd, self.mag, self.err = sigma_clip(self.hjd, self.mag, self.err, sigmas=2.0)
    
    def clip(self, limit = 0.05):
        """
        clip lightcurve at given limit
        """
        m = np.mean(self.mag)
        valid = abs(self.mag-m)<limit
        self.hjd = np.compress(valid, self.hjd)
        self.mag = np.compress(valid, self.mag)
        self.err = np.compress(valid, self.err)
        
    def __len__(self):
        return len(self.hjd)
    
    def psd(self, minperiod, maxperiod):
        from psd import ppsd
        # perform a power spectrum analysis
        tpsa, mpsa = self.hjd, self.mag
        n = len(tpsa)
        # zero padded lightcurves
        t_padded = np.zeros(4*n)
        t_padded[:n] = tpsa
        t_padded[n:] = np.linspace(max(tpsa),4*max(tpsa),3*n)
        m_padded = np.zeros(4*n)
        m_padded[:n] = mpsa
        
        px, f = ppsd(t_padded, 
                     m_padded, 
                     lower=1./maxperiod, 
                     upper=1./minperiod,
                     num= 2000)
        px = np.sqrt(px)
        period = 1./f[np.argmax(px)]
        return period
    
    def pdm(self, minperiod, maxperiod):
        from pdm import pdm
        # look at 20 days or at most at the length of dataset
        import os
        filename = os.path.join(config.lightcurvespath,'pdm',self.starid+'.pdm')
        try:
            
            periods, thetas = np.genfromtxt(filename, unpack=True)
        except IOError:
            periods, thetas = pdm(self.hjd, self.mag, minperiod, maxperiod, 0.5/24)
            a = np.column_stack((periods, thetas))
            np.savetxt(filename, a)
        return periods, thetas
        
    def clean(self, minperiod, maxperiod):
        from clean import clean  # @UnresolvedImport
        import os
        filename = os.path.join(config.lightcurvespath,'clean',self.starid+'.clean')
        try:
            p, cf = np.genfromtxt(filename, unpack=True)
        except IOError:
            t = self.hjd
            x = self.mag
            f, cleaned, _ = clean(t, x, threshold=1e-3)
            n2 = len(f) /2
            cf = cleaned[n2+1:]/(2.0*np.var(x))
            p = 1./f[n2+1:]
            cf = cf[(p>=minperiod) & (p<maxperiod)]
            p = p[(p>=minperiod) & (p<maxperiod)]
            #i = np.argmax(cf)
            #period = p[i]
            a = np.column_stack((p, cf))
            np.savetxt(filename, a)
        return p, cf
            
    def phased(self, period):
        from functions import phase
        tp, yp = phase(self.hjd, self.mag, period)

            
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, resid,_,_ = np.linalg.lstsq(A,yp)
        amp_err = resid[0]

        amp = max(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)-\
              min(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)

        tp1 = np.linspace(0.0, period, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)
        return amp, amp_err
        


class Analysis(object):
    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        self.stars = []
        query = "SELECT starid, bv ,vmag ,vmag_err from ngc6633;"
        self._stars = self.wifsip.query(query)
        
        columns = self.wifsip.columns('ngc6633')
        data_types = self.wifsip.data_types('ngc6633')
        
        for c,d in zip(columns, data_types):
            print c,d
        arraydata = []
        for star in self._stars:
            arraydata.append(tuple(star))
        columns = ['starid', 'bv', 'vmag', 'vmag_err']
        data_types = ['S25', np.float16, np.float16, np.float16]
        self.stars = np.array(arraydata, dtype = zip(columns, data_types))
#        self.age = 10**8.557/1e6 # in Myr from Webda
#        self.ebv = 0.031 # from Webda
#        self.dm = 9.53 # from Webda

    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        query = "SELECT starid, vmag, bv" \
        " FROM ngc6633 " \
        " WHERE vmag<18 " \
        " AND vmag < 10 + bv*5.4 " \
        " AND vmag > 7.5 + bv*5.4" \
       # " AND bv>0.4 " \
        " AND NOT good IS NULL" \
        " ORDER BY vmag LIMIT 100;" 

        result = self.wifsip.query(query)
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]
        self.vmag = [s[1] for s in result]
        self.bv = [s[2] for s in result]
        

    def setperiod(self, starid, period, theta=1.0):
        self.starid = starid
        params = {'starid': starid, 'period': period, 'theta': theta}
        if np.isfinite(period):
            query = """UPDATE ngc6633 SET period = %(period)f, theta = %(theta)f WHERE starid='%(starid)s';""" % params
        else: 
            query = """UPDATE ngc6633 SET period = NULL, theta=NULL WHERE starid='%(starid)s';""" % params
        self.wifsip.execute(query) 

    def setamp(self, starid, amplitude):
        self.starid = starid
        params = {'starid': starid, 'amp': amplitude}
        if np.isfinite(amplitude):
            query = """UPDATE ngc6633 SET amp = %(amp)f WHERE starid='%(starid)s';""" % params
        else: 
            query = """UPDATE ngc6633 SET amp = NULL WHERE starid='%(starid)s';""" % params
        self.wifsip.execute(query) 


    def setbad(self, starid):
        self.starid = starid
        query = """UPDATE ngc6633 SET good = False WHERE starid='%s';""" % starid
        self.wifsip.execute(query) 
        
    #def __setattr__(self, name, value):
    #    params = {'starid': self.starid, 'name': name, 'value': str(value)}
    #    
    #    query = """UPDATE ngc6633 SET %(name)s = %(value)s WHERE starid='%(starid)s';""" % params
    #    self.wifsip.execute(query)
    
    #def __getattribute__(self, name):
    #    params = {'starid': self.starid, 'name': name}
    #    query = "SELECT %(name)s from ngc6633 WHERE starid='%(starid)s';" % params
    #    result = self.wifsip.query(query)
    #    return result[0]

    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        mean = np.mean(self.mag)
        plt.hlines(mean,min(self.hjd),max(self.hjd),linestyle='--')
        plt.xlim(min(self.hjd),max(self.hjd))
        plt.grid()
        plt.errorbar(self.hjd, self.mag, yerr=self.err*0.5, fmt='o')
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])

    def plot_clean(self, periods, amplitudes):
        plt.plot(periods, amplitudes, 'k')
        i = np.argmax(amplitudes)
        period = periods[i]
        plt.axvline(x = period, color='red', alpha=0.5)
        plt.axhline(np.mean(amplitudes), color='b', ls='--')
        plt.axhline(5.*np.mean(amplitudes), color='g', ls='--')
        plt.xlim(min(periods),max(periods))
        plt.minorticks_on()
    
    def plot_pdm(self, periods, thetas):
        #plt.plot(periods, thetas, 'k')
        
        from scipy import signal
        kernel = signal.gaussian(101, 2)
        n = len(thetas)/2
        padded_thetas = np.lib.pad(thetas, n, mode='constant', constant_values=(1.0,1.0))
        smoothed = signal.fftconvolve(padded_thetas, kernel, mode='same')[n:-n]/5
        plt.plot(periods, smoothed, 'k')
        i = np.argmin(smoothed)
        period = periods[i]
        plt.axvline(x = period, color='red', alpha=0.5)
        plt.axhline(0.8, color='b', ls='--')
        
        plt.xlim(min(periods),max(periods))
        plt.ylim(0.0,1.0)
        plt.minorticks_on()
        
    def plot_phase(self, period):    
        from functions import phase
        tp, yp = phase(self.hjd, self.mag, period)

            
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, _,_,_ = np.linalg.lstsq(A,yp)

        tp1 = np.linspace(0.0, 2*period, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)

        plt.scatter(tp, yp, edgecolor='none', alpha=0.75)
        plt.scatter(tp+period, yp, edgecolor='none', alpha=0.75)
        plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                  linestyle='--', linewidth=2)
        plt.xlim(0.0, 2*period)
        plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
        plt.xlabel('P = %.4f' % period)

        
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

        minperiod = 1.3
        maxperiod = 15
        
        for starid,vmag,bv in zip(self.stars,self.vmag,self.bv):
            print '%-24s '% starid,
            try:
                lc = LightCurve(starid)
                
            except IOError:
                logger.error("Can't load lightcurve %s" % starid)
                print 'no lightcurve'
                self.setbad(starid)
                continue
            
            if len(lc)<50:
                logger.warn("%s: not enough datapoints" % starid)
                print 'not enough datapoints'
                continue                    
            # perform a 3sigma clipping
            i = np.where(lc.hjd > 2456762)[0]
            lc.mag = lc.mag[i]
            lc.hjd = lc.hjd[i]
            lc.normalize()
            lc.clip()
            lc.detrend()
            lc.sigma_clip()
            lc.normalize()
            
            self.hjd = lc.hjd
            self.mag = lc.mag
            self.err = lc.err
            
            clean_periods, clean_amplitudes = lc.clean(minperiod, maxperiod)
            pdm_periods, pdm_thetas = lc.pdm(minperiod, maxperiod)
            #period = clean_periods[np.argmax(clean_amplitudes)]

            from scipy import interpolate
            
            i = np.argsort(clean_periods)
            c_periods = clean_periods[i]
            c_amps = clean_amplitudes[i]
            c_periods = np.insert(c_periods, 0 , 0.0)
            c_amps = np.insert(c_amps, 0 , 0.0)
            c_periods = np.insert(c_periods, -1 , maxperiod)
            c_amps = np.insert(c_amps, -1 , 0.0)
            c_int = interpolate.interp1d(c_periods, c_amps)
            
            # use interpolation function returned by `interp1d`
            sum_amp = c_int(pdm_periods)*(1.-pdm_thetas)  
            sum_amp /= max(sum_amp) 
            i = np.argmax(sum_amp)
            period = pdm_periods[i] 
            theta = pdm_thetas[i]
            
            import functions as fx
            
            
            ci = np.argmax(c_amps)
            try:
                pgf_amp, pgf_mean, pgf_sigma = fx.gauss_fit(pdm_periods, 1.-pdm_thetas,  pdm_thetas[i], period, 1.0)
            except:
                pgf_sigma = 0.0
            try:
                cgf_amp, cgf_mean, cgf_sigma = fx.gauss_fit(c_periods, c_amps,  c_amps[ci], c_periods[ci], 1.0)
            except:
                cfg_sigma= 0.0
            
            sigma_limit = 5.3
            theta_limit = 0.8
            try:
                ci = np.argmax(c_amps)
                
                
                params = {'period': period,
                          'period_err': pgf_sigma,
                          'clean_period': c_periods[ci],
                          'clean_amp': c_amps[ci],
                          'clean_sigma': cgf_sigma,
                          'theta': theta,
                          'freq': period,
                          }
                
                keys = ', '.join(params.keys())
                values = ', '.join([str(v) for v in params.values()])
                
                query = """UPDATE ngc6633 SET (%s) = (%s) WHERE starid='%s';""" % (keys, values, starid)
                self.wifsip.execute(query)
            except:
                print 'Cannot store params for starid %s' % starid
                
            mamp = np.mean(clean_amplitudes)
            if max(clean_amplitudes)> sigma_limit*mamp and pdm_thetas[i] < theta_limit:
                print "%.2f %.1f" % (period,   max(clean_amplitudes)/mamp)
                self.setperiod(starid, period)
            else:
                print '< %.1f sigmas or  theta %.2f > %.2f' % (sigma_limit, pdm_thetas[i], theta_limit)
                self.setperiod(starid, np.nan)
                continue           
            
            amp, _ = lc.phased(period)
            self.setamp(starid, amp)
                
            star = {'tab':0, 'bv':bv}

            plt.subplot(411) ##################################################
            plt.title('%s (%d) V = %.2f B-V=%.2f' % (starid, star['tab'], vmag, bv))
            self.plot_lightcurve()
            
            plt.subplot(412) ##################################################
            self.plot_clean(clean_periods, clean_amplitudes)
            

            plt.subplot(413) ##################################################
            plt.plot(pdm_periods, sum_amp,'g')
            plt.axvline(period, color='b')
            self.plot_pdm(pdm_periods, pdm_thetas)
            
            plt.subplot(414) ##################################################
            self.plot_phase(period)
            
#             plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
#             plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
#                      linestyle='--', linewidth=2)
#             plt.xlim(0.0,period)
#             plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
#             plt.xlabel('P = %.4f' % period)
#             plt.grid()
#             #plt.show()
#             period_err = 0.0
#             comment = 'P=%6.3f+-%.3f a=%.3f+-%.4f %.2f' % \
#             (period, period_err, amp, amp_err, theta)
#             plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
#             plt.close()
#                 
#             logger.info( comment)
#             print comment
            if show: plt.show()
            else:
                plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
            plt.close()
        

        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NGC 6633 analysis')
    parser.add_argument('-a', '--analysis', action='store_true', help='create table')
    parser.add_argument('--calibrate', action='store_true', help='create table')
    parser.add_argument('--calibrate2', action='store_true', help='create table')
    parser.add_argument('--create', action='store_true', help='create table')
    parser.add_argument('--clear', action='store_true', help='clear table')
    parser.add_argument('--getframes', action='store_true', help='get frames')
    parser.add_argument('--sigmas', action='store_true', help='plot cmd')
    parser.add_argument('--bv', action='store_true', help='update B-V')
    parser.add_argument('filter', default='V', type=str, help='filter color to process')

    args = parser.parse_args()
    
    if args.analysis: 
        analysis = Analysis()
        analysis.getstars()
        analysis.analysis(show = False)
    if args.calibrate: calibrate('NGC 6633 BVI', args.filter)
    if args.calibrate2:
        from calibrate2 import Calibrate2
        cal = Calibrate2('NGC 6633 rot NW', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot NE', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot SW', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot SE', filtername=args.filter)
        cal.grid()
        
    if args.create or args.clear or args.getframes or args.sigmas or args.bv:
        from photometry import Photometry
        phot = Photometry(objname='NGC 6633 BVI', filtercol=args.filter, dbname='ngc6633')
    if args.create: phot.createtable()
    if args.clear: phot.cleartable()
    if args.getframes: phot.getframes()
    if args.sigmas: phot.update_sigmas()
    if args.bv: phot.update_bv()
