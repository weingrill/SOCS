#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 9, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
import numpy as np

class ImportRef(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def loadfromfile(self, filename=''):
        """
           1-  5  I5    ---     Seq       Sequential identification number
           7-  8  I2    h       RAh       Right ascension (J2000)
          10- 11  I2    min     RAm       Right ascension (J2000)
          13- 18  F6.3  s       RAs       Right ascension (J2000)
          20- 21  I2    deg     DEd       Declination (J2000)
          23- 24  I2    arcmin  DEm       Declination (J2000)
          26- 30  F5.2  arcsec  DEs       Declination (J2000)
          31- 39  F9.3  arcmin  xpos      x-frame coordinate in arcmin
          40- 48  F9.3  arcmin  ypos      y-frame coordinate in arcmin
          49- 57  F9.3  mag     b-y       ?=99.99 b-y color index
          58- 65  F8.4  mag   e_b-y       ?=9.99  error in b-y
          66- 74  F9.3  mag     Vmag      ?=99.99 V magnitude
          75- 82  F8.4  mag   e_Vmag      ?=9.99  error in Vmag
          83- 91  F9.3  mag     m1        ?=99.99 Stroemgren m1 color index
          92- 99  F8.4  mag   e_m1        ?=9.99  error in m1
         100-108  F9.3  mag     c1        ?=99.99 Stroemgren c1 color index
         109-116  F8.4  mag   e_c1        ?=9.99  error in c1
         117-125  F9.3  mag     Hb        ?=99.99 Stroemgren Hbeta color index
         126-133  F8.4  mag   e_Hb        ?=9.99  error in Hbeta
         135-136  A2    ---     Mem       [M NM] Membership label
        """
        names = ['Seq','RAh','RAm','RAs','DEd','DEm','DEs',
                 'xpos','ypos','b-y','e_b-y','Vmag','e_Vmag','m1','e_m1',
                 'c1','e_c1','Hb','e_Hb','Mem']        
        self.data = np.genfromtxt(filename, dtype=None, names=names)
        print self.data[0]
        
ir = ImportRef()
ir.loadfromfile(config.datapath+'table3.dat.gz')