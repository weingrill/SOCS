#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 8, 2013

@author: Joerg Weingrill <jweingrill@aip.de>

Data reduction Class for M48 observation
'''

import config
import logging
from m48star import M48Star            
import pylab as plt
import numpy as np

logging.basicConfig(filename=config.projectpath+'m48_analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')

class M48Analysis(object):
    '''
    Class to implement all methods for analyzing and plotting
    '''
    def __init__(self, path):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        self.stars = []
        self.path = path
        self.age = 10**8.557/1e6 # in Myr from Webda
        self.ebv = 0.031 # from Webda
        self.dm = 9.53 # from Webda
        
    def clearperiods(self):
        """
        reset the periods in the database table
        """
        if not raw_input('press Y to erase the periods in table ')=='Y':
            return
        query="""UPDATE m48stars 
        SET period=NULL,period_err=NULL,theta=NULL,amp=NULL,amp_err=NULL
        WHERE period>0;
        """
        logger.info('resetting periods ...')
        self.wifsip.execute(query)
        
    def getstars(self, allstars=False, maglimit=18.2):
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
#             query = """SELECT starid 
#             FROM m48stars 
#             WHERE bv> 0.4
#             AND vmag<4*bv+13
#             AND vmag < %f
#             AND period IS NULL
#             OR period<0
#             ORDER BY vmag;""" % maglimit
            query = """SELECT starid 
            FROM m48stars 
            WHERE good
            AND vmag < %f
            ORDER BY vmag;""" % maglimit

#             query = """SELECT starid 
#             FROM m48stars 
#             WHERE vmag<4.87*bv+11.6 
#             AND vmag>4.87*bv+10.5;"""
        
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
            f = open(config.resultpath+star+'.tsv', 'wt')
            for s in zip(periods,thetas):
                f.write('%f\t%f\n' % s)
        finally:
            f.close()
    
    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        
        mean = np.mean(self.m)
        #std = np.std(self.m)
        plt.hlines(mean,min(self.t),max(self.t),linestyle='--')
        #plt.ylim(mean+std*3,mean-std*3)
        plt.xlim(min(self.t),max(self.t))
        plt.grid()
        #plt.scatter(self.t, self.m, edgecolor='none')
        plt.errorbar(self.t, self.m, yerr=self.e*0.5, fmt='o')
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])
        
    def analysis(self, show=False):
        """perform a PDM analysis on each lightcurve"""
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

        minperiod = 1.2/24
        maxperiod = 15
        
        for starid in self.stars:
            star = M48Star(starid)
            print '%-24s '% starid,
            try:
                
                t, m, e = star.lightcurve()
                t -= min(t)
            except AttributeError:
                logger.error("Can't load lightcurve %s" % starid)
                print 'no lightcurve'
                continue
            
            if len(t)<50:
                logger.warn("%s: not enough datapoints" % starid)
                print 'not enough datapoints'
                continue                    
            # perform a 3sigma clipping
            self.t, self.m, self.e = sigma_clip(t, m, e)
            
            
            # perform a power spectrum analysis
            tpsa, mpsa = self.t, self.m- np.mean(self.m)
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
            # look at 20 days or at most at the length of dataset
            pdm_periods, pdm_thetas = pdm(self.t, self.m, minperiod, maxperiod, 0.5/24)
            period = pdm_periods[np.argmin(pdm_thetas)]
            psd_period = 1./f[np.argmax(px)]
            psd_freq = f[np.argmax(px)]
            
            from scipy import interpolate
            
            i = np.argsort(1./f)
            psd_periods = 1./f[i]
            psd_power = px[i]
            psd_periods = np.insert(psd_periods, 0, 0.0)
            psd_power = np.insert(psd_power,0 ,0.0)
            psd_int = interpolate.interp1d(psd_periods, psd_power)
            
            # use interpolation function returned by `interp1d`
            sum_amp = psd_int(pdm_periods)*(1.-pdm_thetas)   
            i = np.argmax(sum_amp)
            period = pdm_periods[i] 
            theta = pdm_thetas[i]
            
            star['freq'] = psd_freq
            if abs(star['pman']-star['clean_period'])<star['clean_sigma']:
                star['period'] = star['clean_period']
                star['period_err'] = star['clean_sigma']
            elif abs(star['pman']-period)<np.sqrt(period):
                star['period'] = period
                star['period_err'] = np.sqrt(period)
            else:
                star['period'] = None
                star['period_err'] = None 
            star['theta'] = theta
            
            period = star['period']
            period_err = star['period_err']
            tp, yp = phase(self.t, self.m, period)

                
            s1 = np.sin(2*np.pi*tp/period)
            c1 = np.cos(2*np.pi*tp/period)
            s2 = np.sin(4*np.pi*tp/period)
            c2 = np.cos(4*np.pi*tp/period)
            
            A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
            c, resid,_,_ = np.linalg.lstsq(A,yp)
            star['s1'],star['c1'],star['s2'],star['c2'] = c[1:5]
            amp_err = resid[0]
            star['amp_err'] = amp_err

            amp = max(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)-\
                  min(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)
            star['amp'] = amp

            good = False
            if amp_err<0.09 and amp>0.01 and amp<0.15 and theta<0.72:
                good = True
            else: 
                good = False

            tp1 = np.linspace(0.0, period, 100)
            s1 = np.sin(2*np.pi*tp1/period)
            c1 = np.cos(2*np.pi*tp1/period)
            s2 = np.sin(4*np.pi*tp1/period)
            c2 = np.cos(4*np.pi*tp1/period)
                


            plt.subplot(411) ##################################################
            plt.title('%s (%d) B-V=%.2f' % (starid, star['tab'], star['bv']))
            self.plot_lightcurve()
            
            plt.subplot(412) ##################################################
            plt.axvline(x = psd_period, color='green', alpha=0.5)
            plt.axvline(x = period, color='red', alpha=0.5)
            plt.semilogx(1./f,px*1000, 'k')
            plt.xlim(0.1, 30)
            plt.grid()

            plt.subplot(413) ##################################################
            plt.plot(pdm_periods, pdm_thetas, 'k')
            from functions import normalize
            plt.plot(pdm_periods, normalize(sum_amp), 'b')
            
            #plt.ylim(theta, 1.0)
            plt.axvline(x = psd_period, color='green')
            plt.axvline(x = period, color='red')
            plt.grid()
            
            
            plt.subplot(414) ##################################################
            plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
            plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                     linestyle='--', linewidth=2)
            plt.xlim(0.0,period)
            plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
            plt.xlabel('P = %.4f' % period)
            plt.grid()
            #plt.show()
            comment = 'P=%6.3f+-%.3f a=%.3f+-%.4f %.2f' % \
            (period, period_err, amp, amp_err, theta)
            if show:
                plt.show()
            elif good:
                plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
            plt.close()
                
            logger.info( comment)
            if good: print comment,'*'
            else: print comment
            
    
    def export(self):
        """
        export db to file
        """
        query = """select tab, bv, sqrt(bmag_err^2+vmag_err^2) "bv_err", vmag, 
        vmag_err, clean_period, clean_sigma 
        from m48stars 
        where good 
        order by tab;"""
        comments = """# tab   = table number from publication
# B-V   = B-V color
# B-V_e = error of B-V
# Vmag  = V magnitude
# V_err = error of V magnitude
# P     = period from cleaned spectrum
# P_err = period error width of freq-peak in cleaned spectrum
"""

        header = "#tab  B-V   B-V_e  Vmag   V_err   P      P_err" 

        data = self.wifsip.query(query)
        
        np.savetxt(config.datapath+'periods.txt', 
                   data, 
                   fmt='%4d %.3f %.4f %.3f %.4f %7.3f %.3f',
                   comments = comments,
                   header = header)

        query = """SELECT tab, bv, vmag, member 
        FROM m48stars 
        WHERE NOT bv IS NULL
        ORDER BY tab;"""
        data = self.wifsip.query(query)

        with open(config.datapath+'bvtable.txt','wt') as textfile:
            textfile.write('#tab  B-V    Vmag   mem\n')
            for d in data:
                textfile.write('%4d %6.3f %7.4f %-5.5s\n' % d)

    def export_lightcurves(self):
        for starid in self.stars:
            star = M48Star(starid)
            if star['bv'] > 0.4:
                print '%-24s '% starid,
                try:
                    t, m, e = star.lightcurve()
                    a = np.column_stack((t,m,e))
                except TypeError:
                    print 'no data'
                else:
                    filename = config.lightcurvespath+'%s.dat' % starid
                    np.savetxt(filename, (a), fmt='%.6f %.3f %.4f')
                    print 'exported'
        
     
    def __exit__(self):
        self.wifsip.close()
                
    def tables(self):
        '''
        produce the tables for the publication
        '''
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        
        query = """SELECT tab, vmag, bv, p_fin, e_pfin, amp, member, 
            simbad, notes, provisional
            FROM m48stars 
            WHERE good
            ORDER BY vmag ;"""
        result = self.wifsip.query(query)
        
        from record import Record  # @UnresolvedImport
        datarecord = Record(key='tab', 
            columns = ['tab', 'vmag', 'bv', 'p_fin', 'e_pfin', 'amp', 'member', 
            'simbad', 'notes', 'provisional'])
        datarecord.append(result)
        
        f = open(config.resultpath+'table2.tex','wt')
        f.write("""\\begin{longtable}{rccccccccl}
\caption{\label{tab:rotators}Rotation periods of stars. "p" marks a provisional member and "c" a candidate member}
\hline\hline
Id & V   & B--V & P    & Perr & amp & mem & c/p & BJG \#\\\\
   & mag & mag  & days & days & mag &     &     &     \\\\
\hline
\endfirsthead\n
\caption{continued.}\\
\hline\hline
Id & V   & B--V & P    & Perr  & amp &  mem & c/p & BJG \#\\\\
   & mag & mag  & days & days  & mag &      &     &     \\\\
\hline
\endhead
\hline
\endfoot
%%%%\n""")        
        for d in datarecord:
            if type(d['simbad']) is str and d['simbad'].find('Cl* NGC 2548 ')==0:
                d['simbad'] = d['simbad'][13:]
            if str(d['simbad']) == 'None': d['simbad'] = ''
            if str(d['notes']) == 'None': d['notes'] = ''
            d['memstr'] = '--'
            if d['member']: d['memstr']='M'
            elif d['member']==False: d['memstr']='N'
            
            d['prostr'] = '--'
            if d['provisional']: d['prostr']='p'
            elif not d['provisional']: d['prostr']='c'
             
            try:
                s =  '%(tab)4d & %(vmag).3f & %(bv).2f & %(p_fin).2f & %(e_pfin).2f & %(amp).3f & %(memstr)2s & %(prostr)s & %(simbad)s \\\\ %% %(notes)s\n' % d
                print s,
                f.write(s)
            except TypeError:
                print d
        f.write('\\end{longtable}\n')
        f.write('\\end{longtab}\n')
        f.close()
        return 
        query = """SELECT tab, vmag, vmag_err, bv, bmag_err, ra, dec, member, simbad, notes
            FROM m48stars 
            WHERE not bv IS NULL
            ORDER BY vmag;"""
        data = self.wifsip.query(query)
        f = open(config.resultpath+'table_appendix.tex','wt')
        f.write("""\\begin{longtab}
\\begin{longtable}{rcccccl}
\\caption{\label{tab:appendix}Results of photometric measurements from STELLA WiFSIP.}\\\\
\\hline\\hline
Id &  V  & err & B--V & err & R.A.  & Dec   & mem & simbad \\\\
   & mag & mag & mag  & mag & h:m:s & d:m:s &     & name   \\\\
\hline
\endfirsthead\n
\caption{continued.}\\\
\hline\hline
Id &  V  & err & B--V & err & R.A.  & Dec   & mem & simbad \\\\
   & mag & mag & mag  & mag & h:m:s & d:m:s &     & name   \\\\
\hline
\endhead
\hline
\endfoot
%%%%\n""")
        
        for d in data:
            #print d
            #i = data.index(d)+70
            tab, vmag, vmag_err, bv, bmag_err, ra, dec, member, simbad, notes = d
            
            if type(simbad) is str and simbad.find('Cl* NGC 2548 ')==0:
                simbad = simbad[13:]
            if str(simbad) == 'None': simbad = ''
            if str(notes) == 'None': notes = ''
            memstr='--'
            if member: memstr='M'
            elif member==False: memstr='N'
            try:
                bv_err = np.sqrt(vmag_err**2+bmag_err**2)
            except TypeError:
                bv_err = 0.0 
            c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)  # @UndefinedVariable
            ra_str =  c.ra.to_string(unit=u.hourangle,sep=':', precision=1)  # @UndefinedVariable
            dec_str =  c.dec.to_string(sep=':', precision=0)
            if vmag_err is None: vmag_err=0.0
            if bmag_err is None: bmag_err=0.0
            
            try:
                s =  '%4d & %6.3f & %5.3f & %6.3f & %5.3f & %s & %s & %s & %s \\\\\n' % \
                      (tab, vmag, vmag_err, bv, bv_err, ra_str, dec_str, memstr, simbad)
                if len(notes)>0: 
                    s = s.rstrip('\n')
                    s += ' %% %s\n' % notes
                print s,
                f.write(s)
            except TypeError:
                print d
        f.write('\\end{longtable}\n')
        f.write('\\end{longtab}\n')
        f.close()

    def set_tab_column(self):
        query = """SELECT starid
            FROM m48stars 
            WHERE NOT bv IS NULL
            ORDER BY vmag;"""
        result = self.wifsip.query(query)
        starids = [r[0] for r in result]
        
        for starid in starids:
            tab = starids.index(starid)+1
            print '%4d %s' % (tab,starid)
            query = "UPDATE m48stars set tab=%d WHERE starid='%s';" % (tab,starid)
            self.wifsip.execute(query)
        pass

    def load_periods(self):
        f = open(config.datapath+'periods.txt')
        lines = f.readlines()
        f.close()
        
        for l in lines:
            if l[0]=="#":
                continue
            ls = l.split()
            tab = int(ls[0].strip())
            P_fin = float(ls[10].strip())
            E_pfin = float(ls[11].strip())
            provisional = ls[12] == 'p'
            print tab,P_fin,E_pfin, provisional
            query = """UPDATE m48stars SET (p_fin,e_pfin, provisional)=(%f,%f,%s)
                WHERE tab=%d;""" % (P_fin, E_pfin, provisional, tab)
            self.wifsip.execute(query, commit=False)
        self.wifsip.commit()
        
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='M48 analysis')
    parser.add_argument('--clear', action='store_true', help='clear periods')
    parser.add_argument('-a', '--analysis', action='store_true', help='analysis')
    parser.add_argument('-cal2', action='store_true', help='calibrate lightcurves')
    parser.add_argument('-e', '--export', action='store_true', help='export to textfile')
    parser.add_argument('-l', '--load', action='store_true', help='import textfile')
    parser.add_argument('--allstars', action='store_true', help='fetch all stars')
    parser.add_argument('--lightcurves', action='store_true', help='export lightcurves')
    parser.add_argument('--tables', action='store_true', help='export latex tables')
    parser.add_argument('-tab', action='store_true', help='update tab column')
    
    args = parser.parse_args()
    
    m48 =  M48Analysis(config.datapath)
    if args.clear: m48.clearperiods()
    m48.getstars(allstars=args.allstars, maglimit=18.5)
    
    if args.cal2:
        from calibrate2 import Calibrate2
        fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE','M 48 rot SW']
        for field in fields:
            cal = Calibrate2(field, filtername='V')
            cal.grid()

    #m48.set_simbad()
    if args.analysis: m48.analysis()
    if args.lightcurves: m48.export_lightcurves()
    if args.tab: m48.set_tab_column()
    if args.load: m48.load_periods() 
    if args.tables: m48.tables()
    if args.export: m48.export() 
