#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 13, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
from matplotlib import pyplot as plt


def plot_lightcurve(t, x, p0 = None):
    plt.plot(t, x*1000, 'ko-')
    for p in np.arange(0.0, t[-1], p0):
        plt.axvline(p, color='b', alpha=0.25)
    plt.minorticks_on()
    plt.ylim(max(x*1000)+2, min(x*1000)-2)
    plt.xlim(t[0],t[-1])
    plt.ylabel('mmag')
    
def plot_dft(t, x, p0 = 0.0, axis= None):
    from clean import dft  # @UnresolvedImport
    ft, f = dft(t, x)
    plt.plot(1./f[f>0.0], abs(ft)[f>0.0]/(2.0*np.var(x)), 'k')
    i = np.argmax(ft)
    period = 1./f[i]
    plt.axvline(1, color='r', alpha=0.5, label='1.0 day')
    plt.axvline(p0, color='b', alpha=0.5, label='%.2f days' % p0)
    plt.axvline(period, color='g', alpha=0.5, label='max. peak')

    
    plt.minorticks_on()
    plt.xlim(0.0,max(t)/3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'DFT', va='top', horizontalalignment='right', transform=axis.transAxes)

def plot_window(t, x, axis= None):
    from clean import window  # @UnresolvedImport
    #window
    ftw, fw = window(t, x)
    plt.plot(1./fw[fw>0.0], abs(ftw)[fw>0.0], 'k')
    plt.xlim(0.0,max(t)/2)
    plt.text(0.95, 0.9, 'window', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)

def plot_lomb(t, x, p0 = 0.0, axis= None):
    #lomb scargle
    import scipy.signal as signal
    nout = 100
    flomb = np.linspace(1./60., 1./0.05, nout)
    pgram = signal.lombscargle(t, x, flomb)
    norm = np.sqrt(4*(pgram/nout))/(2.0*np.var(x))
    p = 2.*np.pi/flomb
    plt.plot(p, norm,'k')

    i = np.argmax(norm)
    period = p[i]
    plt.axvline(1, color='r', alpha=0.5, label='1.0 day')
    plt.axvline(p0, color='b', alpha=0.5, label='%.2f days' % p0)
    plt.axvline(period, color='g', alpha=0.5, label='max. peak')

    
    plt.minorticks_on()
    plt.xlim(0.0,max(t)/3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'Lomb-Scargle', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)

def plot_clean(t, x, p0 = 0.0, axis= None):
    from clean import clean  # @UnresolvedImport
    #f, cleaned, _ = clean(t, x, gain=0.9, threshold=2e-3)
    f, cleaned, _ = clean(t, x, threshold=1e-3)
    n2 = len(f) /2
    cf = cleaned[n2+1:]/(2.0*np.var(x))
    p = 1./f[n2+1:]
    cf = cf[(p>=0.5) & (p<23.0)]
    p = p[(p>=0.5) & (p<23.0)]
    i = np.argmax(cf)
    period = p[i]
    
    plt.axvline(1, color='r', alpha=0.5, label='1.0 day')
    plt.axvline(p0, color='b', alpha=0.5, label='%.2f days' % p0)
    plt.axvline(period, color='g', alpha=0.5, label='max. peak')
    plt.xlim(0.0,max(t)/3)
    plt.plot(p, cf, 'k')
    plt.minorticks_on()
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'CLEAN', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)
    
    cf = cf[p>=1.1]
    p = p[p>=1.1]
    i = np.argmax(cf)
    period = p[i]
    return period

def plot_phase(t, x, period):
    from functions import phase
    tp, yp = phase(t, x, period)
    plt.scatter(tp,yp*1000, c='k')
    plt.scatter(tp+period,yp*1000, c='k')
    plt.xlim(0.0,period*2)
    plt.ylim(min(yp)*1000,max(yp)*1000)
    plt.xlabel('period = %.2f' % period)
    plt.axvline(period, linestyle='--', color='k')
    plt.minorticks_on()
    plt.ylabel('mmag')
    plt.ylim(plt.ylim()[::-1])

def plot_pdm(t, x, p0 = 0.0, axis= None):
    from pdm import pdm
    pdm_periods, pdm_thetas = pdm(t, x, 0.1, 30.0, 0.02)
    period = pdm_periods[np.argmin(pdm_thetas)]
    plt.plot(pdm_periods, pdm_thetas, 'k')
    
    plt.axvline(1, color='r', alpha=0.5, label='1.0 day')
    plt.axvline(p0, color='b', alpha=0.5, label='%.2f days' % p0)
    plt.axvline(period, color='g', alpha=0.5, label='max. peak')
    plt.xlim(0.0,max(t)/3)
    plt.minorticks_on()
    plt.ylabel('theta')
    plt.text(0.95, 0.9, 'PDM', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)

    
def do_shuffle(t, x):
    from random import shuffle
    shuffle(x)
    return t, x
    

def detrend(t, x):
    par = np.polyfit(t, x, 1)
    x -= np.polyval(par, t)
    return t, x

def plot_star(tab):
    from m48star import M48Star
    print tab
    
    star = M48Star(None, tab = tab)
    #print star
    lc = star.lightcurve()
    lc.sigma_clip()
    t = lc.hjd
    x = lc.mag
    t -= t[0]
    x -= np.mean(x)
    t, x = detrend(t, x)
    p_fin = star['p_fin']
    print tab, p_fin
    plt.figure()
    plt.subplot('511')
    plt.title('star %d' % tab)
    plot_lightcurve(t, x, p0 = p_fin)
    
    axis = plt.subplot('512')
    #plot_dft(t, x, p0 = p_fin, axis= axis)
    plot_pdm(t, x, p0 = p_fin, axis= axis)
    
    axis = plt.subplot('513')
    plot_lomb(t, x, p0 = p_fin, axis= axis)

    axis = plt.subplot('514')
    period = plot_clean(t, x, p0 = p_fin, axis= axis)

    plt.subplot('515')
    plot_phase(t, x, period)
    
    filename = '/work2/jwe/SOCS/M48/plots/clean_star%d.pdf' % tab 
    print   filename
    plt.savefig(filename)
    plt.close()
    

if __name__ == '__main__':
    """
     284     3.25d     0.08d
     303     1.80d     0.02d
     425     6.28d     0.44d
     517     1.73d     0.01d
     752     8.23d     0.51d
    """
    from matplotlib import rcParams
    params = {'backend': 'Agg',
      'axes.labelsize': 8,
      'axes.titlesize': 10,
      'font.size': 8,
      'xtick.labelsize': 8,
      'ytick.labelsize': 8,
      'figure.figsize': [17/2.54, 24/2.54],
      'savefig.dpi' : 300,
      'font.family': 'sans-serif',
      'axes.linewidth' : 0.5,
      #'xtick.major.size' : 2,
      #'ytick.major.size' : 2,
      }
    rcParams.update(params)
    for tab in [284, 303, 425, 517, 752]:
        plot_star(tab)
    
    
    