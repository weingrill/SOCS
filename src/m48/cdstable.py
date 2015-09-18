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

        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, 
                                 user=config.dbuser, 
                                 host=config.dbhost,
                                 dictcursor=True)
        
    def savetable2(self):
        query = """SELECT tab, vmag, bv, p_fin, e_pfin, amp, amp_err, member, 
            simbad, provisional
            FROM m48stars 
            WHERE good
            ORDER BY vmag ;"""
        self.data = self.wifsip.query(query)
        #Id   Vmag   B-V  P     P_er amp   amp_e M P Simbad
        #         1         2         3         4
        #1234567890123456789012345678901234567890123456789012
        # 287 13.597 0.56  2.88 0.10 0.018 0.004 M p BJG 3157


        
        f = open(config.resultpath+'table2.dat','wt')
        
        
        for d in self.data:
            #i = data.index(d)+1
            simbad = ''
            if type(d['simbad']) is str and d['simbad'].find('Cl* NGC 2548 ')==0:
                simbad = d['simbad'][13:]
            if str('simbad') == 'None': simbad = ''
            memstr = '-'
            if d['member']: memstr='M'
            elif d['member']==False: memstr='N'
            prostr = '-'
            if d['provisional']: prostr='p'
            elif not d['provisional']: prostr='c'
            params = d.copy()
            params['simbad'] = simbad
            params['memstr'] = memstr
            params['prostr'] = prostr 
            s =  '%(tab)4d %(vmag)6.3f %(bv)4.2f %(p_fin)5.2f %(e_pfin).2f %(amp).3f %(amp_err).3f %(memstr)s %(prostr)s %(simbad)s\n' % params
            print s,
            f.write(s)
        f.close()

    def saveappendix(self):
        from numpy import sqrt
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        #TODO: omit simbad column
        query = """SELECT tab, vmag, vmag_err, bv, bmag_err, coord[0] "ra", coord[1] "dec", member, simbad
            FROM m48stars 
            WHERE not bv IS NULL
            ORDER BY vmag;"""
        data = self.wifsip.query(query)
        f = open(config.resultpath+'table_a1.dat','wt')
        
        for d in data:
            simbad = ''
            if type(d['simbad']) is str and d['simbad'].find('Cl* NGC 2548 ')==0:
                simbad = d['simbad'][13:]
            simbad = simbad.rstrip()
            memstr='-'
            if d['member']: memstr='M'
            elif d['member']==False: memstr='N'
            try:
                bv_err = sqrt(d['vmag_err']**2+d['bmag_err']**2)
            except TypeError:
                bv_err = 0.0 
            c = SkyCoord(ra=d['ra']*u.deg, dec=d['dec']*u.deg, frame='icrs')  # @UndefinedVariable
            ra_str =  c.ra.to_string(unit=u.hourangle, sep=' ', pad= True, precision=2)  # @UndefinedVariable
            dec_str =  c.dec.to_string(sep=' ', precision=2, pad= True)
            posstr = c.to_string(precision=5)
            #if bmag_err is None: bmag_err=0.0
            params = d.copy()
            params['memstr'] = memstr
            params['bv_err'] = bv_err
            params['ra_str'] = ra_str
            params['dec_str'] = dec_str
            params['posstr'] = posstr
            if d['vmag_err'] is None: params['vmag_err'] = 0.0
            if str(d['simbad']) == 'None': params['simbad'] = ''
            
            try:
                s =  '%(tab)4d %(posstr)s %(vmag)6.3f %(vmag_err)5.3f %(bv)6.3f %(bv_err)5.3f %(ra_str)s %(dec_str)s %(memstr)s %(simbad)-15s \n' % params
                print s,
                f.write(s)
            except TypeError:
                print d, len(s)
        f.close()

    def update_coords(self, filtercol = 'V'):
        import numpy as np
        from tools import log
        query = """SELECT starid, ra, dec 
        FROM m48stars 
        WHERE NOT bv IS NULL AND coord IS NULL ORDER BY tab;"""
        result = self.wifsip.query(query)
        stars = []
        for starid, ra, dec in result:
            stars.append({'starid': starid, 'ra': ra, 'dec': dec})
            
        for star in stars:
            #print star['NUMBER'],star['MAG_ISO'], 
            params = {'ra': star['ra'],
                      'dec': star['dec'],
                      'filtercol': filtercol}
            query = """SELECT alphawin_j2000, deltawin_j2000 
                FROM phot, frames
                WHERE object like 'M 48 %%'
                AND phot.objid = frames.objid
                AND filter = '%(filtercol)s'
                AND circle(phot.coord,0) <@ circle(point(%(ra).11f,%(dec).11f), 0.3/3600.0)
                AND flags<5;""" % params
            result = self.wifsip.query(query)
            #print result
            if len(result) == 0:
                print '%-25s: no data' % star['starid']
                query = """UPDATE m48stars 
                SET coord = NULL
                WHERE starid = '%s';""" % star['starid']
                self.wifsip.execute(query)
                continue
            elif len(result)==1:
                print '%-25s: not enough data' % star['starid']
                query = """UPDATE m48stars 
                SET coord = NULL
                WHERE starid = '%s';""" % star['starid']
                self.wifsip.execute(query)
                continue
            elif len(result)>1:
                ra = np.array([r[0] for r in result])
                dec = np.array([r[1] for r in result])
                
                mean_ra = np.mean(ra)
                mean_dec = np.mean(dec)
                std_ra = np.std(ra)
                std_dec = np.std(dec)
                # remove outliers
                std = np.sqrt(std_ra**2 + std_dec**2)
                dist = np.sqrt((ra - mean_ra)**2 + (dec - mean_dec)**2) 
                i = np.where(dist < 2*std)
                mean_ra = np.mean(ra[i])
                mean_dec = np.mean(dec[i])
                std_ra = np.std(ra[i])
                std_dec = np.std(dec[i])
                n = len(ra[i])
            params['starid'] = star['starid']
            params['ra'] = mean_ra
            params['dec'] = mean_dec
            params['n'] = n
            
            params['sra'] = std_ra*3600.0
            params['sdec'] = std_dec*3600.0
            query = """UPDATE m48stars 
            SET coord = point(%(ra).11f,%(dec).11f)
            WHERE starid = '%(starid)s';""" % params
            #query = query.replace('nan', 'NULL')
            self.wifsip.execute(query)
            log(config.projectpath+'coords1.log','%(starid)-25s, %(ra).6f, %(dec).6f, %(sra).2f, %(sdec).2f, %(n)3d' % params)
        self.wifsip.commit()    
        

if __name__ == '__main__':
    cdst = CDSTable() 
    #cdst.update_coords()
    #cdst.savetable2()
    cdst.saveappendix()
