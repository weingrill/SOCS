#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 14, 2014

@author: jwe
'''

from m48star import M48Star
import pylab as plt
from functions import sigma_clip
import numpy as np

starids = ['20140301A-0071-0013#125',
           '20140303A-0074-0013#1855',
           '20140301A-0071-0013#1218',
           '20140301A-0071-0013#1316',
           '20140306A-0097-0013#732',
           '20140307A-0000-0013#334']

plt.subplots(3, 2)
for starid in starids:            
    star = M48Star(starid)
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
    
    ax = plt.subplot(3,2,starids.index(starid))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.hlines(mean,min(t),max(t),linestyle='--')
    plt.xlim(min(t),max(t))
    plt.plot(t,m,'gray')
    plt.scatter(t, m, edgecolor='none', facecolor='k', s=20)
    ylim=[max(m)+0.01, min(m)-0.01]
    plt.text(1, ylim[1]+0.0025, star['tab'], fontsize=12)
    plt.ylim(ylim[1],ylim[0])
    #plt.xlabel('days')
    #plt.ylabel('V mag')
plt.tight_layout()    
#plt.show()
plt.savefig('/work2/jwe/m48/plots/lightcurves.pdf')
plt.savefig('/work2/jwe/m48/plots/lightcurves.png')
    