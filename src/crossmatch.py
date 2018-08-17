#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 29, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
from matplotlib import pyplot as plt
from datasource import DataSource

from astropy.coordinates import SkyCoord, search_around_sky  # @UnresolvedImport
from astropy import units as u


class Frame():

    def __init__(self, objid):
        self.objid = objid
        table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
        query = """SELECT * FROM frames WHERE objid = '%s';""" % objid
        result = table.query(query)[0]
        print(result.keys())
        self.exposure = result['expt']


def getcoords(objid):
    table = DataSource(database='stella', user='stella', host='pera.aip.de')
    query = """SELECT alphawin_j2000, deltawin_j2000, flux_auto 
                FROM phot
                WHERE objid = '%s'
                ORDER BY star;""" % objid
    result = table.query(query)
    ra = np.array([r[0] for r in result])
    dec = np.array([r[1] for r in result])
    flux = np.array([r[2] for r in result])
    return ra, dec, flux


ra1, dec1, f1 = getcoords('20130602A-0003-0016')  # b
ra2, dec2, f2 = getcoords('20130602A-0003-0018')  # y

frame1 = Frame('20130602A-0003-0016')
frame2 = Frame('20130602A-0003-0018')
f1 /= frame1.exposure
f2 /= frame2.exposure

# get imaging data
# image_data = fetch_imaging_sample()
c = SkyCoord(ra=ra1 * u.degree, dec=dec1 * u.degree)
catalog = SkyCoord(ra=ra2 * u.degree, dec=dec2 * u.degree)
# idx, d2d, d3d = c.match_to_catalog_sky(catalog)
# matches = catalog[idx]
# print len(matches)

# idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1.5*u.arcsec)
idxc, idxcatalog, d2d, d3d = search_around_sky(c, catalog, 1.5 * u.arcsec)
print(len(idxc))
# plt.hist(d2d*3600.0, bins=30)
# plt.show()
mby = -2.5 * np.log10(f1[idxc] / f2[idxcatalog])
my = -2.5 * np.log10(f2[idxcatalog])
plt.scatter(mby, my)
plt.ylim(plt.ylim()[1], plt.ylim()[0])
plt.grid()
plt.show()

plt.close()
