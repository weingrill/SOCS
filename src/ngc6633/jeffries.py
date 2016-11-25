#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 17, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np

def avg(jid, rvs, errs):
    rvs = np.array(rvs)
    errs = np.array(errs)
    wghts = 1./errs
    #print 'J%d %.1f %.1f' % (jid, np.average(rvs, weights=wghts), np.linalg.norm(errs) )
    
    return np.average(rvs, weights=wghts), np.linalg.norm(errs) 

"""    
avg(65, [-28.2, -34.6], [2.4, 2.3])
avg(66, [-27.8,-32.6,-30.4],[2.3,2.3,2.0])
avg(69, [-25.7,35.0],[5.1,4.7])
avg(70, [-26.7,-27.5],[4.4,4.3])
avg(71, [-34.4, -35.0],[2.2,2.0])
avg(74, [-27.3, -29.5], [3.2, 2.6])
avg(82, [-32.1, -25.5, -26.4], [2.9, 2.3, 2.3])
avg(83, [-31.2, -25.9, -26.7], [2.2, 2.0, 2.0])
avg(88, [-32.5, -35.4], [4.2, 3.3])
avg(92, [-29.7, -28.1], [2.2, 2.1])
avg(103, [-23.7, -30.0], [2.2, 2.0])
"""

with open('/work2/jwe/SOCS/NGC6633/data/jeffries2002.txt') as f:
    lines = f.readlines()

from astropy import units as u
from astropy.coordinates import SkyCoord  # @UnresolvedImport


jeffriesdb = {}
for l in lines[1:]:
    ls = l.split()
    #Ident rh rm rs dd dm ds V B-V V-I RV RVerr
    jid = int(ls[0].lstrip('J'))
    rh,rm,rs,dd,dm,ds,v,bv,vi = ls[1:10]
    
    
    starcoordinates =  SkyCoord('%s %s %s %s %s %s' %(rh,rm,rs,dd,dm,ds), unit=(u.hourangle, u.deg))    
    
    rv = float(ls[10])
    rverr = float(ls[11])
    if not jid in jeffriesdb:
        jeffriesdb[jid] = {'ra': starcoordinates.ra.deg,
                           'dec': starcoordinates.dec.deg,
                           'V': float(v),
                           'B-V': float(bv),
                           'V-I': float(vi),
                           'rv':[rv], 'rverr': [rverr]}
    else:
        jeffriesdb[jid]['rv'].append(rv)
        jeffriesdb[jid]['rverr'].append(rverr)

for j in jeffriesdb:
    rv, rverr = avg(j, jeffriesdb[j]['rv'], jeffriesdb[j]['rverr'])
    #print j, rv, rverr
print jeffriesdb