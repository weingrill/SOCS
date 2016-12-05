#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 30, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import config
from datasource import DataSource
wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
#import matplotlib.pyplot as plt
import numpy as np

def fprint(s):
    with open(config.datapath+'/ngc6633phot.log', 'at') as f:
        f.write(s)

def calibrate(objid, filtercol):
    ffield = {'V': '"V"', 'B': '"B"', 'I': '"I"'}
    queryparams = {'objid': objid, 
                   'ffield': ffield[filtercol]}
    query = """SELECT mag_isocor, magerr_isocor, %(ffield)s FROM phot, referencestars 
                WHERE objid = '%(objid)s'
                AND circle(phot.coord,0) <@ circle(referencestars.coord, 1.0/3600.0)
                AND flags < 4 AND mag_isocor>0 AND %(ffield)s > 0""" % queryparams
    stars = wifsip.query(query)
    if len(stars) < 15:
        print '0 0 0 0 0 0'
        fprint('0 0 0 0 0 0\n')
        return
    mag_isocor = np.array([r[0] for r in stars])
    #magerr_isocor = np.array([r[1] for r in stars])
    refmag = np.array([r[2] for r in stars])
    ocmag = mag_isocor - refmag
    std = np.std(ocmag)
    print '%d %f %f' % (len(stars), np.median(ocmag), std),
    fprint('%d %f %f ' % (len(stars), np.median(ocmag), std))
    mag_isocor1 = mag_isocor - np.median(ocmag)
    ocmag1 = mag_isocor1 - refmag
    k = abs(ocmag1 - std) < std
    if len(k) < 3:
        print '0 0 0'
        fprint('0 0 0\n')
        return
    mag_isocor2 = mag_isocor[k]
    refmag2 = refmag[k]
    mean =  np.mean(mag_isocor2 - refmag2)
    
    print '%d %f %f' % (len(mag_isocor2), mean, np.std(mag_isocor2 - refmag2))
    fprint('%d %f %f\n' % (len(mag_isocor2), mean, np.std(mag_isocor2 - refmag2)))
    #plt.title(objid)
    #plt.plot(refmag2, mag_isocor2-refmag2-mean, 'ob')
    #plt.xlabel('%s' % ffield[filtercol])
    #plt.show()

fprint('#objid filter expt field stars median std stars2 mean std2\n')
for filtercol in ['V', 'B']:
    for expt in [60, 300]:
        if filtercol == 'B': expt *=2
        for field in ['C', 'NW', 'NE', 'SW', 'SE']:
            queryparams = {'filtercol': filtercol, 
                           'expt': str(expt),
                           'field': field}
            query = "SELECT objid FROM frames " + \
                    " WHERE object LIKE 'NGC 6633 BVI %(field)s' AND filter='%(filtercol)s' AND expt=%(expt)s;" % queryparams
            #print query
            objids = wifsip.query(query)
            for objid in objids:
                print objid[0], filtercol, expt, field,
                fprint('%s %s %d %s ' % (objid[0], filtercol, expt, field))
                calibrate(objid[0], filtercol)
            
            