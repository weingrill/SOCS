#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 4, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

    
from datasource import DataSource
from glob import glob
import config
import os

lightcurves = glob(config.lightcurvespath+'20140625A-0014-0001*.dat')

wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
def findstar(starid):
    objid, star = starid.split('#')
    query = """SELECT alphawin_j2000, deltawin_j2000
    FROM phot WHERE objid='%s' AND star=%d; """ % (objid, int(star))
    result = wifsip.query(query)
    ra, dec = result[0]
    query = """SELECT starid
    FROM ngc6633
    WHERE circle(point(ra, dec),0) <@ circle(point(%.11f, %.11f),1.0/3600); """ % (ra, dec)
    result = wifsip.query(query)
    #print starid, ra, dec, result
    if result:
        return result[0][0]
    else:
        return starid
 
for lc in lightcurves:
    lcname = lc[lc.rfind('/')+1:]
    starid = lcname.rstrip('.dat')
    
    newname = config.lightcurvespath+findstar(starid)+'.dat'
    print lc, newname
    if lc<>newname:
        os.rename(lc, newname)
wifsip.close()