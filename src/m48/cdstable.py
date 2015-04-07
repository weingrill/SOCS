#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 1, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config

class CDSTable(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.datastruct = []
    def savetable2(self):
        from datasource import DataSource
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT tab, vmag, bv, p_fin, e_pfin, amp, amp_err, member, 
            simbad, provisional
            FROM m48stars 
            WHERE good
            ORDER BY vmag ;"""
        self.data = wifsip.query(query)
        d0 =  self.data[0]
        print d0
        #Id   Vmag   B-V  P     P_er amp   amp_e M P Simbad
        #         1         2         3         4
        #1234567890123456789012345678901234567890123456789012
        # 287 13.597 0.56  2.88 0.10 0.018 0.004 M p BJG 3157


        
        f = open(config.resultpath+'table2.dat','wt')
        
        
        for d in self.data:
            #print d
            #i = data.index(d)+1
            tab, vmag, bv, period, period_err, \
            amp, amp_err,member, simbad, provisional = d

            if type(simbad) is str and simbad.find('Cl* NGC 2548 ')==0:
                simbad = simbad[13:]
            if str(simbad) == 'None': simbad = ''
            memstr = '-'
            if member: memstr='M'
            elif member==False: memstr='N'
            prostr = '-'
            if provisional: prostr='p'
            elif not provisional: prostr='c'
             
            s =  '%4d %6.3f %4.2f %5.2f %.2f %.3f %.3f %s %s %s\n' % \
            (tab, vmag, bv, period, period_err, amp, amp_err, memstr, prostr, simbad)
            print s
            f.write(s)
        f.close()

    def saveappendix(self):
        from numpy import sqrt
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        from datasource import DataSource
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        
        query = """SELECT tab, vmag, vmag_err, bv, bmag_err, ra, dec, member, simbad
            FROM m48stars 
            WHERE not bv IS NULL
            ORDER BY vmag;"""
        data = wifsip.query(query)
        f = open(config.resultpath+'table_appendix.dat','wt')
        for d in data:
            #print d
            #i = data.index(d)+70
            tab, vmag, vmag_err, bv, bmag_err, ra, dec, member, simbad = d
            
            if type(simbad) is str and simbad.find('Cl* NGC 2548 ')==0:
                simbad = simbad[13:]
            if str(simbad) == 'None': simbad = ''
            memstr='-'
            if member: memstr='M'
            elif member==False: memstr='N'
            try:
                bv_err = sqrt(vmag_err**2+bmag_err**2)
            except TypeError:
                bv_err = 0.0 
            c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)  # @UndefinedVariable
            ra_str =  c.ra.to_string(unit=u.hourangle,sep=' ', precision=1)  # @UndefinedVariable
            dec_str =  c.dec.to_string(sep=' ', precision=0)
            if vmag_err is None: vmag_err=0.0
            if bmag_err is None: bmag_err=0.0
            
            try:
                s =  '%4d %6.3f %5.3f %6.3f %5.3f %s %s %s %s \n' % \
                      (tab, vmag, vmag_err, bv, bv_err, ra_str, dec_str, memstr, simbad)
                print s,
                f.write(s)
            except TypeError:
                print d
        f.close()

if __name__ == '__main__':
    cdst = CDSTable() 
    cdst.savetable2()
    cdst.saveappendix()
