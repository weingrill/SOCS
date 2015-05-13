#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 18, 2014

@author: jwe

class definiton for paper plots on M48
'''
from m48star import M48Star            
import numpy as np
import pylab as plt
from matplotlib import rcParams
from functions import sigma_clip, phase

class M48Plots(object):
    '''
    class definiton for paper plots on M48
    '''

    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='oldpina.aip.de')
        self.stars = []
        self.getstars()

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

    def plot_lightcurve(self, starid):
        """
        plot the lightcurve for a given star
        """
        star = M48Star(starid)
        try:
            t, m, e = star.lightcurve()
        except TypeError:
            print 'no data'
            return
        t, m, e = sigma_clip(t, m, e)
        t -= t[0]
        mean = np.mean(m)
        m -= mean
        mean = 0.0
        plt.hlines(mean,min(t),max(t),linestyle='--')
        plt.xlim(min(t),max(t))
#        plt.grid()
        plt.scatter(t, m, edgecolor='none', facecolor='k', s=5)
        plt.plot(t,m,'gray')
        #plt.errorbar(t, m, yerr=e*0.5, fmt='o', s=5)
        #ylim=plt.ylim()
        ylim=[max(m)+0.01, min(m)-0.01]
        plt.text(1, ylim[1]+0.0025, star['tab'], fontsize=12)
        plt.ylim(ylim[1],ylim[0])

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
            ax = plt.subplot(7,4,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self.plot_lightcurve(starid)
            if sp==7*4+1 or starid==self.stars[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig('/work2/jwe/m48/plots/lightcurves%d.pdf' % lc)
                    plt.savefig('/work2/jwe/m48/plots/lightcurves%d.eps' % lc)
                lc += 1
                sp = 1 
                plt.close()


    def phase_plot(self, starid, axis):
        star = M48Star(starid)
        try:
            t, m, e = star.lightcurve()
        except TypeError:
            print 'no data'
            return
        t, m, e = sigma_clip(t, m, e)
        m -= np.mean(m)
        period = star['p_fin']
        tp, yp = phase(t, m, period)
        plt.axvline(period, linestyle='-.', color='0.5')
        plt.axhline(0.0, linestyle='--', color='0.5')
        plt.scatter(tp, yp-np.mean(yp), edgecolor='none', facecolor='k', s=5)
        plt.scatter(tp+period, yp-np.mean(yp), edgecolor='none', facecolor='k', s=5)
        plt.xticks(np.arange(60))
        plt.xlim(0,period*2)
        plt.ylim(plt.ylim()[1],plt.ylim()[0])
        plt.text(0, 0, star['tab'], 
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
            ax = plt.subplot(7,4,sp)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            sp += 1
            self.phase_plot(starid, ax)
            if sp==7*4+1 or starid==self.stars[-1]:
                plt.tight_layout()
                if show: plt.show()
                else: 
                    plt.savefig('/work2/jwe/m48/plots/phase%d.pdf' % phase)
                    plt.savefig('/work2/jwe/m48/plots/phase%d.eps' % phase)
                phase += 1
                sp = 1 
                plt.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='M48 analysis')
    parser.add_argument('-p','--phaseplots', action='store_true', 
                        help='make phaseplots')
    parser.add_argument('-l','--lightcurves', action='store_true', 
                        help='make lightcurve plots')
    parser.add_argument('-show', action='store_true', 
                        help='show plots instead of saving to file')

    args = parser.parse_args()
    
    m48plots = M48Plots()
    if args.phaseplots: m48plots.make_phaseplots(show=args.show)
    if args.lightcurves: m48plots.make_lightcurves(show=args.show)