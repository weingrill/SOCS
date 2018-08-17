#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 17, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from opencluster import OpenCluster
from astropy.coordinates import SkyCoord  # @UnresolvedImport

if __name__ == '__main__':
    iriscoord = SkyCoord('02h05m4.4s', '+21d24m06.1s', frame='icrs')

    Iris = OpenCluster(objectname='Iris',
                       obsmode='uvby',
                       ra=iriscoord.ra.degree,
                       dec=iriscoord.dec.degree)
    Iris.startdate = '2017-10-30'
    Iris.enddate = '2017-10-31'
    Iris.tofile('/work2/jwe/SOCS/')
