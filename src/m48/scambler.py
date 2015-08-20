#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 18, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from matplotlib import pyplot as plt
#from psd import psd, ppsd
from clean import clean  # @UnresolvedImport
import numpy as np
from random import shuffle

def _init(argst, argsy):
    global t_global
    global y_global
        
    t_global = argst
    y_global = argsy

def _worker(wi):
    global t_global
    global y_global
    
    shuffle(y_global)
    #px, f = psd(t_global, y_global, lower=1./15., upper=1./0.05)
    f, px, _, sigma0 = clean(t_global, y_global)
    px = px[f>=0.0]
    f = f[f>=0.0]
    i = np.argmax(px)
    if f[i] ==  0.0:
        return np.inf, 0.0
    else:
        return 1./f[i], px[i]/sigma0


def test(doshuffle = True):
    from m48star import M48Star
    
    star = M48Star(None, tab=284)
    lc = star.lightcurve()
    t = lc.hjd
    t-=t[0]
    y = lc.mag
    y -= np.mean(y)
    if doshuffle: 
        shuffle(y)
    f, px, _, sigma0 = clean(t, y)
    px = px[f>=0.0]
    f = f[f>=0.0]
    
    i = np.argmax(px)
    print '%.2f %.1f' % (1./f[i], px[i]/sigma0)
    
    plt.plot(1./f[1:], px[1:]/sigma0, 'r')
    plt.xlim(0.0,15.0)
    plt.axhline(5, linestyle='--')
    plt.show()

def pooltest():
    from m48star import M48Star
    
    star = M48Star(None, tab=284)
    lc = star.lightcurve()
    t  = lc.hjd
    t -= t[0]
    y  = lc.mag
    y -= np.mean(y)
    
    par = np.polyfit(t, y, 1)
    
    y -= par[0]*t + par[1]

    f, px, _, sigma0 = clean(t, y)
    px = px[f>=0.0]
    f = f[f>=0.0]
    
    i = np.argmax(px)
    fi = f[i]
    pxi = px[i]
    print '%.2f %.1f' % (1./fi, pxi/sigma0)
    
    runs = 10000
    
    from multiprocessing import Pool
    w = np.arange(runs)
    pool = Pool(initializer=_init, initargs=(t,y))
    p = pool.map(_worker, w)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks
    f = np.array([pi[0] for pi in p])
    px = np.array([pi[1] for pi in p])
    plt.subplot('211')
    k = np.argsort(f)
    f = f[k] 
    px = px[k]
    np.savetxt('/work2/jwe/SOCS/M48/data/scrambler.txt', np.c_[f,px], fmt='%6.3f',header='period sigma')
    px1 = px[(f>3.07) & (f<3.33)]
    print px1.shape, 1./fi, pxi
    plt.scatter(f, px, edgecolor='none', alpha=0.5)
    plt.scatter(1./fi, pxi/sigma0, c='r', edgecolor='none', s=80)
    plt.minorticks_on()
    #plt.show()
    plt.xlim(0.0, t[-1]/2)
    plt.title('star 284: %d runs' % runs)
    plt.xlabel('period [days]')
    plt.ylabel('$\sigma$')
    plt.subplot('212')
    plt.hist(f, bins=np.sqrt(runs), range=[0.0,t[-1]/2], normed=True)
    plt.axvline(1./fi, color='r')
    plt.savefig('/work2/jwe/SOCS/M48/plots/scrambler1.pdf')
    
    plt.close()


if __name__ == '__main__':
    #test(doshuffle=False)
    pooltest()