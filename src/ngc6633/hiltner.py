#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 24, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

with open('/work2/jwe/SOCS/NGC6633/data/hiltner.txt') as f:
    lines = f.readlines()

print 'lines read: %d' % len(lines)

#from astropy import units as u
#from astropy.coordinates import SkyCoord  # @UnresolvedImport


hiltnerarray = []
for i,l in enumerate(lines[1:]):
    ls = l.split()
    #Number x y V B-V U-B n Member
    #print i, ls
    if len(ls)>8 or len(ls)<7:
        print 'wrong number of elements in line %d: "%s"' % (i+1, l.rstrip('\n'))
        
    try:
        number = int(ls[0])
    except ValueError:
        print 'invalid value %s for Number in line %d' % (ls[0], i+1)
    try:
        x = float(ls[1])
    except ValueError:
        print 'invalid value %s for x in line %d' % (ls[1], i+1)
    try:
        y = float(ls[2])
    except ValueError:
        print 'invalid value %s for y in line %d' % (ls[2], i+1)
    try:
        v = float(ls[3])
    except ValueError:
        print 'invalid value %s for V in line %d' % (ls[3], i+1)
    try:
        bv = float(ls[4])
    except ValueError:
        print 'invalid value %s for B-V in line %d' % (ls[4], i+1)
    try:
        ub = float(ls[5])
    except ValueError:
        print 'invalid value %s for U-B in line %d' % (ls[5], i+1)
    try: 
        n = int(ls[6])
    except ValueError:
        print 'invalid value %s for n in line %d' % (ls[6], i+1)
    if len(ls)==8:
        member = str(ls[7])
        if member not in ['X','D','X?']:
            print 'invalid value %s for member in line %d' % (ls[7], i+1)
    else: member=''
    
    hiltnerarray.append((number, x, y, v, bv, ub, n, 0., 0.))

import numpy as np
    
hilt = np.array(hiltnerarray, dtype = [('number', int), ('x', float), ('y', float), ('V', 'float'), ('B-V', 'float'), ('U-B','float'), ('n', int),('R.A.', float), ('Dec', float)])

refarray = []
from astropy import units as u
from astropy.coordinates import SkyCoord  # @UnresolvedImport

# taken from http://simbad.harvard.edu/simbad/sim-ref?simbo=on&submit=submit%20bibcode&querymethod=bib&bibcode=1958ApJ...127..539H
for n, cstr in [(28, '18 26 47 +06 30.3'),
                (35, '18 26 54 +06 32.6'),
                (47, '18 27 02 +06 30.8'),
                (49, '18 27 03 +06 27.4'),
                (53, '18 27 04 +06 37.7'),
                (54, '18 27 05 +06 34.6'),
                (76, '18 27 22 +06 27.3'),
                (80, '18 27 24 +06 27.9'),
                (95, '18 27 38 +06 37.3'),
                (99, '18 27 39 +06 37.7'),
                (101, '18 27 43 +06 37.1'),
                (149, '18 28 43.015 +06 48 47.52'),
                (151, '18 29 04.1303 +06 47 48.660'),
                (154, '18 29 14.6408 +06 52 08.036'),
                (158, '18 29 23.425 +06 35 06.24'),
                (161, '18 29 58.323 +06 46 28.67')]:    
    c =  SkyCoord(cstr, unit=(u.hourangle, u.deg))
    refarray.append((n,c.ra.deg, c.dec.deg))  

ref = np.array(refarray, dtype = [('number', int), ('R.A.', float), ('Dec', float)])

ns = ref['number']
xs = hilt['x'][ns-1]
ys = hilt['y'][ns-1]
ra = ref['R.A.']
dec = ref['Dec']
A = np.vstack([xs, ys, np.ones(len(ns))]).T
wcsra = np.linalg.lstsq(A, ra)[0]
wcsdec = np.linalg.lstsq(A, dec)[0]
print wcsra
print wcsdec
dx1, dy1, cx = wcsra
dx2, dy2, cy = wcsdec
hilt['R.A.'] = dx1*hilt['x'] + dy1*hilt['y'] + cx
hilt['Dec'] = dx2*hilt['x'] + dy2*hilt['y'] + cy
np.savetxt('/work2/jwe/SOCS/NGC6633/data/hiltner_solved.txt', hilt, fmt='%03d %5.2f %5.2f %5.2f %-6.2f %-6.2f %d %.6f %6f', header='#Number x y V B-V U-B n ra dec')

n, xs, ys, ra, dec = np.loadtxt('/work2/jwe/SOCS/NGC6633/data/hiltner_apass.csv', delimiter=',', unpack=True)
A = np.vstack([xs, ys, np.ones(len(n))]).T
wcsra = np.linalg.lstsq(A, ra)[0]
wcsdec = np.linalg.lstsq(A, dec)[0]
print wcsra
print wcsdec
hilt['R.A.'] = dx1*hilt['x'] + dy1*hilt['y'] + cx
hilt['Dec'] = dx2*hilt['x'] + dy2*hilt['y'] + cy
np.savetxt('/work2/jwe/SOCS/NGC6633/data/hiltner_solved2.txt', hilt, fmt='%03d %5.2f %5.2f %5.2f %-6.2f %-6.2f %d %.6f %6f', header='#Number x y V B-V U-B n ra dec')


