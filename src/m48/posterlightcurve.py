#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 14, 2014

@author: jwe
'''

if __name__ == '__main__':
    from m48star import M48Star
    import pylab as plt
    from functions import sigma_clip
    import numpy as np            
    star = M48Star('20140303A-0074-0013#1855')
    try:
        t, m, e = star.lightcurve()
    except TypeError:
        print 'no data'
        exit()
    t, m, e = sigma_clip(t, m, e)
    t -= t[0]
    mean = np.mean(m)
    m -= mean
    mean = 0.0
    plt.hlines(mean,min(t),max(t),linestyle='--')
    plt.xlim(min(t),max(t))
    plt.plot(t,m,'gray')
    plt.scatter(t, m, edgecolor='none', facecolor='k', s=20)
    ylim=[max(m)+0.01, min(m)-0.01]
    #plt.text(1, ylim[1]+0.0025, star['tab'], fontsize=12)
    plt.ylim(ylim[1],ylim[0])
    plt.xlabel('days')
    plt.ylabel('V mag')
    #plt.show()
    plt.savefig('/work1/jwe/Dropbox/Documents/Poster/AIP2014/m48_lc332.eps')
    