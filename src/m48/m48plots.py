#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 18, 2014

@author: jwe

class definiton for paper plots on M48
'''
import config
from m48star import M48Star            
import numpy as np
import pylab as plt
from matplotlib import rcParams

class M48Plots(object):
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
        self.age = 10**8.557/1e6 # in Myr from Webda
        self.ebv = 0.031 # from Webda
        self.dm = 9.53 # from Webda

    def getstars(self):
        """
        build up a list of stars, that have been marked sa good
        """
        
        query = """SELECT starid 
        FROM m48stars 
        WHERE good
        ORDER BY tab;"""
        
        result = self.wifsip.query(query)
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]

    def _plot_lightcurve(self, starid, axis):
        """
        plot the lightcurve for a given star
        """
        star = M48Star(starid)
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
        plt.yticks(np.arange(-0.1,0.1,0.01))
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
        star = M48Star(starid)
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
        
        star = M48Star(starid)
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
        import astronomy as ast
        from dbtable import DBTable  # @UnresolvedImport
        
        astars = DBTable(self.wifsip, 'm48stars', condition='vmag<12.0 and good is NULL')
        gstars = DBTable(self.wifsip, 'm48stars', condition='good and provisional is null')
        pstars = DBTable(self.wifsip, 'm48stars', condition='provisional')
        
        fig_width = 8.9/2.54  # width in inches, was 7.48in
        fig_height = 8.9/2.54  # height in inches, was 25.5
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
          'axes.linewidth' : 0.2,
          'xtick.major.size' : -2,
          'ytick.major.size' : -2,
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
        ax.scatter(pstars['ra'],pstars['dec'], marker='o',edgecolor='r', facecolor='none', s=10)
        ax.scatter(gstars['ra'],gstars['dec'], marker='o',edgecolor='r', facecolor='r', s=10)
        
        ra = np.array([8.24242, 8.21705, 8.20416])*15.0
        de = np.array([-6.08887,-5.70876,-5.51204])
        
        x = [ra[0],ra[2],ra[2],ra[1],ra[1],ra[0],ra[0]]
        y = [de[0],de[0],de[1],de[1],de[2],de[2],de[0]]
        
        
        cra, cdec = (123.42916666666666,-5.75)
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
        xticks = [ast.hms2dd((8,m,0)) for m in [16,15,14,13,12]]
        xlabels = ['$8^h%2d^m$' % m for m in [16,15,14,13,12]]
        plt.xticks(xticks, xlabels)
        declist = [(-5,15),(-5,30),(-5,45),(-6,0),(-6,15)]
        yticks = [ast.dms2dd(dl) for dl in declist]
        ylabels = ['$%d^{\circ}%2d^m$' % dl for dl in declist]
        plt.yticks(yticks, ylabels, rotation=-90)
        ax.grid()
        ax.set_ylim(cdec-0.4,cdec+0.4)
        ax.set_xlim(cra+0.4,cra-0.4)
        plt.xlabel('R.A. (J2000)')
        plt.ylabel('Dec. (J2000)')
        if show:
            plt.show()  
        else:
            plt.savefig(config.plotpath+'m48_map.pdf', transparent=True)
            plt.savefig(config.plotpath+'m48_map.eps', transparent=True)

    def plot_cmd(self, show=False, mark_active=False):
        from dbtable import DBTable  # @UnresolvedImport
        
        alls = DBTable(self.wifsip, 'm48stars', condition='NOT bv is NULL')
        good = DBTable(self.wifsip, 'm48stars', condition='good and provisional is null')
        prov = DBTable(self.wifsip, 'm48stars', condition='provisional')
        
        iso_v, iso_bv = self.load_isochrone() 
        plt.plot(iso_bv, iso_v, 'g', alpha=0.3, lw=5.0,label='800 Myr iso')

        plt.scatter(alls['bv']-self.ebv,alls['vmag'], edgecolor='none', alpha=0.75, s=4, c='k')

        if mark_active:
            plt.scatter(good['bv']-self.ebv, good['vmag'], marker='o',edgecolor='r', facecolor='r', s=30, label='rotators')
            plt.scatter(prov['bv']-self.ebv, prov['vmag'], marker='o',edgecolor='r', facecolor='none', s=30, label='provisional')
        plt.legend()
        plt.title('M48 Color Magnitude Diagram')
        plt.ylim(20.0, 10.0)
        plt.xlim(0.0, 2.0)
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('V [mag]')
        plt.grid()
        if show:
            plt.show()
        else:
            if mark_active:
                plt.savefig(config.plotpath+'m48cmd_active.eps')
                plt.savefig(config.plotpath+'m48cmd_active.pdf')
                plt.savefig(config.plotpath+'m48cmd_active.png', dpi=300)
            else:
                plt.savefig(config.plotpath+'m48cmd.eps')
                plt.savefig(config.plotpath+'m48cmd.pdf')
        plt.close()

    def plot_cpd(self, show=False):
        query = """SELECT bv, period, amp 
                    FROM m48stars 
                    WHERE NOT good
                    AND NOT bv IS NULL
                    AND period>0;"""
        data = self.wifsip.query(query)
        
        bv = np.array([d[0] for d in data])
        period = np.array([d[1] for d in data])
        
        query = """SELECT bv, period, amp 
                    FROM m48stars 
                    WHERE vmag<4*bv + 13 
                    AND NOT bv IS NULL
                    AND period>0
                    AND good;"""
        data = self.wifsip.query(query)
        bv_ms = np.array([d[0] for d in data])
        period_ms = np.array([d[1] for d in data])
        amp_ms = np.array([d[2] for d in data])
        
        query = """SELECT bv, period, amp 
                    FROM m48stars 
                    WHERE NOT bv IS NULL
                    AND period>0
                    AND member;"""
        data = self.wifsip.query(query)
        bv_mem = np.array([d[0] for d in data])
        period_mem = np.array([d[1] for d in data])
        amp_mem = np.array([d[2] for d in data])
        
        
        import gyroage
        from functions import logspace
        
        bvgyro = logspace(0.5, 2.0, num=100)
        #P = gyroage.gyroperiod(bv360, 360.0, version=2007)
        P, pc = gyroage.gyroperiod(bvgyro, 360.0, version=2010)
        plt.plot(bvgyro, pc, color='b', linestyle='--')
        plt.plot(bvgyro, P, color='r')
        
        plt.scatter(bv-self.ebv, period, s=1, 
                    edgecolor='none', c='k')
        plt.scatter(bv_ms-self.ebv, period_ms, s=(1.0-amp_ms)*50., 
                    edgecolor='none', facecolor='green')
        
        plt.scatter(bv_mem-self.ebv, period_mem, s=(1.0-amp_mem)*50., 
                    edgecolor='blue',
                    facecolor='none')
        
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('period [days]')
        plt.ylim(0.0, 20.0)
        plt.xlim(0.0, 2.0)
        plt.grid()
        if show:
            plt.show()
        else:
            plt.savefig(config.plotpath+'m48cpd.eps')
            plt.savefig(config.plotpath+'m48cpd.pdf')
        plt.close()

    def load_isochrone(self):
        from numpy import loadtxt
        isofile = config.datapath+'0p500Gyr_FeH0p0_Y0p277_AMLTsol.iso'
        a = loadtxt(isofile)
        iso_mv = a[:,5]
        iso_bv = a[:,6]
        return iso_mv+self.dm, iso_bv
    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='M48 analysis')
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
    
    m48plots = M48Plots()
    if args.phaseplots: m48plots.make_phaseplots(show=args.show)
    if args.lightcurves: m48plots.make_lightcurves(show=args.show)
    if args.spectra: m48plots.make_spectra(show=args.show)
    if args.map: m48plots.plot_map(show=args.show)
    if args.cmd: 
        m48plots.plot_cmd(show=args.show)
        m48plots.plot_cmd(show=args.show, mark_active=True)
    if args.cpd: m48plots.plot_cpd(show=args.show)
    