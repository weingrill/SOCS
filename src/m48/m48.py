#!/usr/bin/python
'''
Created on May 8, 2013

@author: jwe <jweingrill@aip.de>

Data reduction Class for M48 observation
'''

import logging
logging.basicConfig(filename='/work2/jwe/m48/m48_analysis.log', 
                    format='%(asctime)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')



from m48star import M48Star            

class M48Analysis(object):
    '''
    classdocs
    '''


    def __init__(self, path):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.stars = []
        self.path = path
        self.age = 10**8.557/1e6 # in Myr from Webda
        self.ebv = 0.031 # from Webda
        self.dm = 9.53 # from Webda
    
    def clearperiods(self):
        """
        reset the periods in the database table
        """
        
        query="""UPDATE m48stars 
        SET period=NULL,period_err=NULL,theta=NULL,amp=NULL,amp_err=NULL
        WHERE period>0;
        """
        logger.info('resetting periods ...')
        self.wifsip.execute(query)
        
    def getstars(self, allstars=False, maglimit=21.0):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        if allstars:
            query = """SELECT starid 
            FROM m48stars 
            WHERE NOT bv IS NULL
            AND period IS NULL
            ORDER BY vmag;"""
        else:
            query = """SELECT starid 
            FROM m48stars 
            WHERE NOT bv IS NULL
            AND vmag<4*bv+13
            AND vmag < %f
            ORDER BY vmag;""" % maglimit
        
        logger.info('fetching stars ...')
        result = self.wifsip.query(query)
        logger.info('... %d stars found' % len(result))
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]
    
    def set_simbad(self):
        from PySimbad import simcoo
        for star in self.stars:
            ra, dec = self.wifsip.query("""SELECT ra,dec 
            from m48stars 
            where starid='%s'""" % star)[0]
            print ra,dec,
            simbad = simcoo(ra, dec)
            print simbad
            if simbad=='None':
                self.wifsip.execute("""UPDATE m48stars 
                SET simbad=NULL 
                WHERE starid='%s'""" % (simbad, star))
            else:
                self.wifsip.execute("""UPDATE m48stars 
                SET simbad='%s' 
                WHERE starid='%s'""" % (simbad, star))
                
    
    
    def store_pdm(self, star, periods, thetas):
        """
        store the periods and thetas for a given star in a tab separated file
        """
        try:
            f = open('/work2/jwe/m48/results/'+star+'.tsv', 'wt')
            for s in zip(periods,thetas):
                f.write('%f\t%f\n' % s)
        finally:
            f.close()
    
    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        import pylab as plt
        import numpy as np
        
        mean = np.mean(self.m)
        std = np.std(self.m)
        plt.hlines(mean,min(self.t),max(self.t),linestyle='--')
        #plt.ylim(mean+std*3,mean-std*3)
        plt.xlim(min(self.t),max(self.t))
        plt.grid()
        #plt.scatter(self.t, self.m, edgecolor='none')
        plt.errorbar(self.t, self.m, yerr=self.e*0.5, fmt='o')
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])
        
                    
    def analysis(self):
        """perform a PDM analysis on each lightcurve"""
        import numpy as np
        import pylab as plt
        from pdm import pdm
        from psd import ppsd
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
            star = M48Star(starid)
            print '%-24s '% starid,
            try:
                
                t, m, e = star.lightcurve()
                t -= min(t)
                
                # perform a 3sigma clipping
                self.t, self.m, self.e = sigma_clip(t, m, e)
                
                
                # perform a power spectrum analysis
                px, f = ppsd(self.t, 
                             self.m- np.mean(self.m), 
                             lower=1./30, 
                             upper=1./0.1,
                             num= 2000)
                px = np.sqrt(px)
                # look at 20 days or at most at the length of dataset
                length = min([max(self.t), 20.0])
                p1, t1 = pdm(self.t, self.m, 0.1, length, 60./86400.)
                period = p1[np.argmin(t1)]
                period1 = 1./f[np.argmax(px)]
                freq = f[np.argmax(px)]
                if abs(period * freq - 1) < 0.25:
                    isgood = True
                else:
                    isgood = False
                star.freq = freq
                theta = min(t1)
                star.period = period
                star.theta = theta
                star.period_err = abs(period-1./freq) 
                amp = np.max(px)
                    
                if np.std(px)>0 and amp>0:
                    star.amp = amp
                    amp_err = amp/np.std(px)
                    star.amp_err = amp_err
                    
                tp, yp = phase(self.t,self.m, period)
            except (TypeError):
                comment = 'no data'
                star.period = None
                star.theta = None
            else:
                plt.subplot(411)
                plt.title('%s B-V=%.2f' % (starid, star.bv))
                self.plot_lightcurve()
                
                plt.subplot(412)
                plt.axvline(x = period1, color='green', alpha=0.5)
                plt.axvline(x = period, color='red', alpha=0.5)
                plt.semilogx(1./f,px, 'k')
                plt.xlim(0.1, 30)
                plt.grid()
    
                plt.subplot(413)
                plt.plot(p1, t1, 'k')
                plt.ylim(theta, 1.0)
                plt.axvline(x = period1, color='green')
                plt.axvline(x = period, color='red')
                plt.grid()
                
                s1 = np.sin(2*np.pi*tp/period)
                c1 = np.cos(2*np.pi*tp/period)
                s2 = np.sin(4*np.pi*tp/period)
                c2 = np.cos(4*np.pi*tp/period)
                
                A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
                c, resid,rank,sigma = np.linalg.lstsq(A,yp)
                star.s1,star.c1,star.s2,star.c2 = c[1:5]
                amp_err = resid[0]/np.sqrt(len(self.m))
                star.amp_err
                plt.subplot(414)
                plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
                tp1 = np.linspace(0.0, period, 100)
                s1 = np.sin(2*np.pi*tp1/period)
                c1 = np.cos(2*np.pi*tp1/period)
                s2 = np.sin(4*np.pi*tp1/period)
                c2 = np.cos(4*np.pi*tp1/period)
                
                plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                         linestyle='--', linewidth=2)
                plt.xlim(0.0,period)
                plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
                plt.grid()
                #plt.show()
                comment = '%6.3f %6.3f %4.2f %.3f %.4f' % \
                (period, period1, theta, amp, amp_err)

                plt.savefig('/work2/jwe/m48/plots/%s.pdf' % starid)
                plt.close()
                
            logger.info( comment)
            print comment
            
    def make_cmd(self, show=False, mark_active=False):
        import pylab as plt
        import numpy as np
        query = "SELECT vmag+0.06, bv FROM m48stars WHERE NOT bv is NULL;"
        data = self.wifsip.query(query)
        vmag = np.array([d[0] for d in data])
        bv = np.array([d[1] for d in data])
        
        iso_v, iso_bv = self.load_isochrone() 
        query = """SELECT vmag+0.06, bv 
                    FROM m48stars 
                    WHERE period>0
                    AND abs(period*freq-1)<0.5 
                    ;"""
        data = self.wifsip.query(query)
        vmag_good = np.array([d[0] for d in data])
        bv_good = np.array([d[1] for d in data])
        print self.ebv
        plt.scatter(bv-self.ebv,vmag, edgecolor='none', alpha=0.75, s=4, c='k')
        plt.plot(iso_bv, iso_v, 'r')
        
        if mark_active:
            plt.scatter(bv_good-self.ebv,vmag_good, edgecolor='none', alpha=0.75, s=30, c='g')
        
        #k = 4
        #d = 13
        #x = np.linspace(-0.5, 2.5, 10)
        #y = k*x+d
        #plt.plot(x, y, linestyle='dashed', color='b')
        plt.ylim(21.0, 8.0)
        plt.xlim(-0.2, 2.0)
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('V [mag]')
        plt.grid()
        if show:
            plt.show()
        else:
            if mark_active:
                plt.savefig('/work2/jwe/m48/plots/m48cmd_active.eps')
                plt.savefig('/work2/jwe/m48/plots/m48cmd_active.pdf')
            else:
                plt.savefig('/work2/jwe/m48/plots/m48cmd.eps')
                plt.savefig('/work2/jwe/m48/plots/m48cmd.pdf')
        plt.close()

    def make_cpd(self, show=False):
        import pylab as plt
        import numpy as np
        query = """SELECT bv, period, theta 
                    FROM m48stars 
                    WHERE NOT good
                    AND NOT bv IS NULL
                    AND period>0
                    AND abs(period*freq-1)<0.5;"""
        data = self.wifsip.query(query)
        
        bv = np.array([d[0] for d in data])
        period = np.array([d[1] for d in data])
        
        query = """SELECT bv, period, theta 
                    FROM m48stars 
                    WHERE vmag<4*bv + 13 
                    AND NOT bv IS NULL
                    AND period>0
                    AND good;"""
        data = self.wifsip.query(query)
        bv_ms = np.array([d[0] for d in data])
        period_ms = np.array([d[1] for d in data])
        theta_ms = np.array([d[2] for d in data])
        
        import gyroage
        from functions import logspace
        
        bv360 = logspace(0.5, 2.0, num=100)
        #P = gyroage.gyroperiod(bv360, 360.0, version=2007)
        P, pc = gyroage.gyroperiod(bv360, 360.0, version=2003)
        plt.plot(bv360, pc, color='b', linestyle='--')
        plt.plot(bv360, P, color='r')
        
        plt.scatter(bv-self.ebv, period, s=1, 
                    edgecolor='none', c='k')
        plt.scatter(bv_ms-self.ebv, period_ms, s=(1.0-theta_ms)*50., 
                    edgecolor='none',
                    alpha=0.75, facecolor='green')
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('period [days]')
        plt.ylim(0.0, 20.0)
        plt.xlim(0.0, 2.0)
        plt.grid()
        if show:
            plt.show()
        else:
            plt.savefig('/work2/jwe/m48/plots/m48cpd.eps')
            plt.savefig('/work2/jwe/m48/plots/m48cpd.pdf')
        plt.close()
    
    def export(self):
        """
        export db to file
        """
        import numpy as np
        query = """SELECT bv, vmag+0.06, period, freq
        FROM m48stars
        WHERE NOT bv IS NULL AND period>0.0"""
        data = self.wifsip.query(query)
        np.savetxt('/work1/jwe/Dropbox/M48/data/periods.txt', 
                   data, 
                   fmt='%.3f %6.3f %.3f %.6f')
        
    def load_isochrone(self):
        from numpy import loadtxt
        isofile = '/work2/jwe/m48/data/0p500Gyr_FeH0p0_Y0p277_AMLTsol.iso'
        a = loadtxt(isofile)
        iso_mv = a[:,5]
        iso_bv = a[:,6]
        return iso_mv+self.dm, iso_bv
     
    def __exit__(self):
        self.wifsip.close()
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='M48 analysis')
    parser.add_argument('--clear', action='store_true', help='clear periods')
    parser.add_argument('-a', '--analysis', action='store_true', help='analysis')
    parser.add_argument('-cmd', action='store_true', help='plot cmd')
    parser.add_argument('-cpd', action='store_true', help='plot cpd')
    parser.add_argument('-e', '--export', action='store_true', help='export to textfile')
    parser.add_argument('--allstars', action='store_true', help='fetch all stars')
    
    args = parser.parse_args()
    
    m48 =  M48Analysis('/work2/jwe/m48/data/')
    if args.clear: m48.clearperiods()
    m48.getstars(allstars=args.allstars, maglimit=21)
    
    #m48.set_simbad()
    if args.analysis: m48.analysis()
    if args.cmd: 
        m48.make_cmd()
        m48.make_cmd(mark_active=True)
    if args.cpd: m48.make_cpd()
    if args.export: m48.export() 
