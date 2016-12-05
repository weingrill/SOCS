#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 18, 2014

@author: jwe

class definiton for paper plots on NGC6633 derived from M48
'''
import config
from ngc6633star import NGC6633Star            
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

class NGC6633Plots(object):
    '''
    class definiton for paper plots on M48
    '''

    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        self.stars = []
        self.getstars()
        self.rows = 7
        self.columns = 3
        self.age = 600 # from Jeffries 2002
        self.ebv = 0.165 # from Jeffries 2002
        self.dm = 7.77 # from Jeffries 2002

    def getstars(self):
        """
        build up a list of stars, that have been marked sa good
        """
        from dbtable import DBTable
        query = """SELECT starid 
        FROM ngc6633 
        WHERE good
        ORDER BY tab;"""
        
        result = self.wifsip.query(query)
        print '... %d stars found' % len(result)
        self.stars = DBTable(self.wifsip, 'ngc6633', condition='good')

    def _plot_lightcurve(self, starid, axis):
        """
        plot the lightcurve for a given star
        """
        star = NGC6633Star(starid)
        try:
            lc = star.lightcurve()
            lc.sigma_clip()
            lc.rebin(0.01)
            t = lc.hjd
            m = lc.mag
        except TypeError:
            print 'no data'
            return
        t -= t[0]
        mean = np.mean(m)
        m -= mean
        mean = 0.0
        plt.axhline(mean,linestyle='--', color='b')
        plt.xlim(min(t),max(t))
        plt.xticks(np.arange(0.0,70.0,10.0))
        plt.yticks(np.arange(-0.2,0.2,0.01))
        plt.scatter(t, m, edgecolor='none', facecolor='k', s=5)
        plt.plot(t,m,'gray')
        ylim=[max(m)+0.01, min(m)-0.01]
        if star['provisional']:
            startab = '(%d)' % star['tab']
        else:
            startab = star['tab']
        plt.text(0.01, 0.01, startab, 
                 fontsize=12, 
                 horizontalalignment='left',
                 verticalalignment='bottom',
                 transform=axis.transAxes)
        plt.ylim(ylim)

    def make_lightcurves(self, show=False):
        """plot lightcurve"""
        
        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)
        
        sp = 1
        lc = 1
        for starid in self.stars:
            print starid
            ax = plt.subplot(self.rows, self.columns,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self._plot_lightcurve(starid, ax)
            if sp==self.rows*self.columns + 1 or starid==self.stars[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig(config.plotpath+'lightcurves%d.pdf' % lc)
                    plt.savefig(config.plotpath+'lightcurves%d.eps' % lc)
                lc += 1
                sp = 1 
                plt.close()


    def phase_plot(self, starid, axis):
        from functions import phase
        star = NGC6633Star(starid)
        try:
            lc = star.lightcurve()
            lc.sigma_clip()
            lc.rebin(0.01)
            t = lc.hjd
            m = lc.mag
        except TypeError:
            print 'no data'
            return
        m -= np.mean(m)
        period = star['p_fin']
        tp, yp = phase(t, m, period)
        plt.axvline(period, linestyle='-.', color='0.5')
        plt.axhline(0.0, linestyle='--', color='0.5')
        plt.scatter(tp, yp-np.mean(yp), edgecolor='none', facecolor='k', s=5)
        plt.scatter(tp+period, yp-np.mean(yp), edgecolor='none', facecolor='k', s=5)
        
        plt.xticks(np.arange(60))
        plt.xlim(0, period*2)
        plt.yticks(np.arange(plt.ylim()[0],plt.ylim()[1],0.01))
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        
        if star['provisional']:
            startab = '(%d)' % star['tab']
        else:
            startab = star['tab']
        
        plt.text(0.01, 0.01, startab, 
                 fontsize=12,
                 verticalalignment='bottom',
                 transform=axis.transAxes)

    def make_phaseplots(self, show=False):
        """plot lightcurve"""
        
        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)

        
        sp = 1
        phase=1
        for starid in self.stars:
            print starid
            ax = plt.subplot(self.rows, self.columns,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self.phase_plot(starid, ax)
            if sp==self.rows * self.columns+1 or starid==self.stars[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig(config.plotpath+'phase%d.pdf' % phase)
                    plt.savefig(config.plotpath+'phase%d.eps' % phase)
                phase += 1
                sp = 1 
                plt.close()

    def _plot_spectrum(self, starid, axis):
        from numpy import std
        
        star = NGC6633Star(starid)
        try:
            freq, amp = star.cleanspectrum()
        except TypeError:
            print 'no data'
            return
        period = star['p_fin']
        
        plt.axvline(period, linestyle='-.', color='green')
        plt.axvline(1.0, linestyle='-.', color='red')
        
        plt.axhline(5.0*std(amp), linestyle='--', color='blue')
        
        plt.plot(1./freq,amp, 'k')
        
        plt.xticks(np.arange(60))
        plt.xlim(0,15.0)
        #plt.ylim(plt.ylim()[1],plt.ylim()[0])
        
        if star['provisional']:
            startab = '(%d)' % star['tab']
        else:
            startab = star['tab']
        plt.text(0.99, 0.99, startab, 
                 fontsize=12,
                 verticalalignment='top',
                 horizontalalignment='right',
                 transform=axis.transAxes)


    def make_spectra(self, show=False):
        """
        plots the CLEANed spectra
        """
        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)
        
        sp = 1
        lc = 1
        for starid in self.stars:
            print starid
            ax = plt.subplot(self.rows, self.columns,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self._plot_spectrum(starid, ax)
            if sp==self.rows*self.columns+1 or starid==self.stars[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig(config.plotpath+'spectra%d.pdf' % lc)
                    plt.savefig(config.plotpath+'spectra%d.eps' % lc)
                lc += 1
                sp = 1 
                plt.close()

    def plot_map(self, show=False):
        '''
        plots the map of the M48 observation including the BJG field
        '''
        #import pywcsgrid2
        import astro.coordinates as ast
        from dbtable import DBTable  # @UnresolvedImport
        
        astars = DBTable(self.wifsip, 'ngc6633', condition='vmag<12.0 and good is NULL')
        gstars = DBTable(self.wifsip, 'ngc6633', condition='good')
        #pstars = DBTable(self.wifsip, 'NGC6633Star', condition='provisional')
        
        #set plot attributes
        plt.style.use('aanda.mplstyle')
        params = {
          u'figure.subplot.bottom' : 0.11,                                                                                                                             
          u'figure.subplot.left' : 0.11,                                                                                                                               
          u'figure.subplot.right' : 0.89,                                                                                                                              
          u'figure.subplot.top' : 0.89
          }
        rcParams.update(params)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        vmag = np.array(astars['vmag'])
        
        ax.set_aspect(1.)
        ax.scatter(astars['ra'],astars['dec'], marker='*', s=(12.-vmag)*15,facecolor='gray', edgecolor='none')
        #ax.scatter(pstars['ra'],pstars['dec'], marker='o',edgecolor='r', facecolor='none', s=10)
        ax.scatter(gstars['ra'],gstars['dec'], marker='o',edgecolor='r', facecolor='r', s=10)
        
        ra = np.array([8.24242, 8.21705, 8.20416])*15.0
        de = np.array([-6.08887,-5.70876,-5.51204])
        
        x = [ra[0],ra[2],ra[2],ra[1],ra[1],ra[0],ra[0]]
        y = [de[0],de[0],de[1],de[1],de[2],de[2],de[0]]
        
        
        cra, cdec = (276.88, 6.57)
        d2 = 0.5*1320.2/3600.0
        
        fields = [(cra    , cdec),
                 (cra - d2, cdec + d2),
                 (cra + d2, cdec + d2),
                 (cra - d2, cdec - d2),
                 (cra + d2, cdec - d2)]
        for field in fields:
            ra,dec = field
            
            ras = [ra-d2, ra+d2, ra+d2, ra-d2, ra-d2]
            das = [dec-d2, dec-d2, dec+d2, dec+d2, dec-d2]
            ax.plot(ras, das ,'g')

        ax.plot(x,y, 'k', linestyle='--',label='BJG 2005')
        
        #xticks = ax.get_xticks()
        print [ast.dd2hms(r) for r in ras]
        xticks = [ast.hms2dd((18,m,0)) for m in [30,29,28,27,26,15]]
        xlabels = ['$8^h%2d^m$' % m for m in [30,29,28,27,26,15]]
        plt.xticks(xticks, xlabels)
        declist = [(6,0),(6,15),(6,30),(6,45),(7,0)]
        yticks = [ast.dms2dd(dl) for dl in declist]
        ylabels = ['$%d^{\circ}%2d^m$' % dl for dl in declist]
        plt.yticks(yticks, ylabels, rotation=90)
        ax.grid()
        ax.set_ylim(cdec-0.4,cdec+0.4)
        ax.set_xlim(cra+0.4,cra-0.4)
        plt.xlabel('R.A. (J2000)')
        plt.ylabel('Dec. (J2000)')
        if show:
            plt.show()  
        else:
            plt.savefig(config.plotpath+'ngc6633_map.pdf', transparent=True)
            plt.savefig(config.plotpath+'ngc6633_map.eps', transparent=True)

    def plot_cmd(self, show=False, mark_active=False, isochrone=True):
        from dbtable import DBTable
        vc = [ 0.98976977,  0.01785897,  0.07519376]
        bc = [ 0.99106353,  0.01346018,  0.07998822]

        alls = DBTable(self.wifsip, 'ngc6633', condition='NOT bv is NULL AND vmag<20')
        good = DBTable(self.wifsip, 'ngc6633', condition='good')
        members = DBTable(self.wifsip, 'ngc6633', condition='member')
        #prov = DBTable(self.wifsip, 'ngc6633', condition='provisional')
        
        vmag = vc[0]*alls['vmag'] + vc[1]*alls['bv'] + vc[2]
        bmag = bc[0]*alls['bmag'] + bc[1]*alls['bv'] + bc[2]
        bv = bmag - vmag
        
        plt.style.use('aanda.mplstyle')
        suffix = ''
        if isochrone:
            print 'Av = %f' % (3.1*self.ebv)
            iso_v, iso_bv = self.load_isochrone() 
            plt.plot(iso_bv+self.ebv, iso_v+self.dm+(3.1*self.ebv), 'b', alpha=0.3, lw=5.0,label='600 Myr iso') #
            suffix = suffix + '_iso'

        plt.scatter(bv, vmag, edgecolor='none', alpha=0.75, s=2, c='k')
        
        if mark_active:
            plt.scatter(good['bv'], good['vmag'], marker='o',edgecolor='r', facecolor='r', s=20, label='rotators')
            suffix = suffix + '_active'
        
        plt.scatter(members['bv'], members['vmag'], marker='o', edgecolor='g', facecolor='none', s=20, label='members')
            
        plt.legend(loc='lower left')
        plt.title('NGC6633 Color Magnitude Diagram')
        plt.ylim(20.0, 7.0)
        plt.xlim(0.0, 1.7)
        plt.xlabel('B - V')
        plt.ylabel('V [mag]')
        plt.minorticks_on()
        plt.grid()
        if show:
            plt.show()
        else:
            plt.savefig(config.plotpath+'ngc6633cmd%s.eps' % suffix)
            plt.savefig(config.plotpath+'ngc6633cmd%s.pdf' % suffix)
        plt.close()

    def plot_cpd(self, show=False, gyrochrone=False):
        from dbtable import DBTable
        
        members = DBTable(self.wifsip, 'ngc6633', condition='member')
        good = DBTable(self.wifsip, 'ngc6633', condition='good and period>0')
        nonmembers = DBTable(self.wifsip, 'ngc6633', condition='NOT member')
        
        
        import gyroage
        from functions import logspace
        
        bvgyro = logspace(0.5, 1.5, num=100)
        age = 550
        P = gyroage.gyroperiod(bvgyro, age, version=2010)
        
        plt.style.use('aanda.mplstyle')
        
        if gyrochrone:
            plt.plot(bvgyro, P, color='b', label='%d Myr'% age)
        
        plt.scatter(good['bv']-self.ebv, good['period'], marker='o',edgecolor='r', facecolor='r', s=20, label='rotators')
        plt.scatter(members['bv']-self.ebv, members['period'], marker='o',edgecolor='k', facecolor='r', s=20, label='members')
        plt.scatter(nonmembers['bv']-self.ebv, nonmembers['period'], marker='x',edgecolor='k', facecolor='k', s=20, label='nonmembers')
        
        
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('period [days]')
        plt.ylim(0.0, 15.0)
        plt.xlim(0.4, 1.6)
        plt.legend()
        plt.title('NGC 6633 Color Period Diagram')
        #plt.grid()
        if show:
            plt.show()
        else:
            if gyrochrone:
                plt.savefig(config.plotpath+'ngc6633cpd_gyro.eps')
                plt.savefig(config.plotpath+'ngc6633cpd_gyro.pdf')
            else:
                plt.savefig(config.plotpath+'ngc6633cpd.eps')
                plt.savefig(config.plotpath+'ngc6633cpd.pdf')
        plt.close()

    def load_isochrone(self):
        from numpy import loadtxt
        isofile = config.datapath+'yapsi_w_X0p70322_Z0p01877_600Myr.dat'
        mass, logT, logL, logg, Mv, ub, bv, vr, vi = loadtxt(isofile, unpack=True)
        return Mv, bv
    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='NGC6633 analysis')
    parser.add_argument('-p','--phaseplots', action='store_true', 
                        help='make phaseplots')
    parser.add_argument('-l','--lightcurves', action='store_true', 
                        help='make lightcurve plots')
    parser.add_argument('-s','--spectra', action='store_true', 
                        help='make lightcurve plots')
    parser.add_argument('-m','--map', action='store_true', 
                        help='plot the map')
    parser.add_argument('-cmd', action='store_true', help='plot cmd')
    parser.add_argument('-cpd', action='store_true', help='plot cpd')
    parser.add_argument('--show', action='store_true', 
                        help='show plots instead of saving to file')

    args = parser.parse_args()
    
    ngc6633plots = NGC6633Plots()
    if args.phaseplots: ngc6633plots.make_phaseplots(show=args.show)
    if args.lightcurves: ngc6633plots.make_lightcurves(show=args.show)
    if args.spectra: ngc6633plots.make_spectra(show=args.show)
    if args.map: ngc6633plots.plot_map(show=args.show)
    if args.cmd: 
        ngc6633plots.plot_cmd(show=args.show)
        ngc6633plots.plot_cmd(show=args.show, mark_active=True)
    if args.cpd: 
        ngc6633plots.plot_cpd(show=args.show)
        ngc6633plots.plot_cpd(show=args.show, gyrochrone=True)
    
