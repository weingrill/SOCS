#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 3, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
import numpy as np

class NGC6633cat(object):
    '''
    fill the ngc6633table in the database
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
    
    def load_positions(self, filename=config.datapath+'NGC6633V.fit'):
        import pyfits
        
        hdulist = pyfits.open(filename)
        data = hdulist[2].data
        i = np.argsort(data['MAG_ISO'])
        self.stars = data[i]
        hdulist.close()
        print self.stars.columns
        print "%d stars loaded" % len(self.stars)
         
        
    def build(self, filtercol = 'V'):
        from psycopg2 import IntegrityError
        for star in self.stars: 
            #print star['NUMBER'],star['MAG_ISO'], 
            params = {'ra': star['ALPHAWIN_J2000'],
                      'dec': star['DELTAWIN_J2000'],
                      'filtercol': filtercol}
            query = """SELECT phot.objid || '#' || phot.star, mag_isocor, corr, magerr_isocor 
                FROM phot, frames
                WHERE object like 'NGC 6633 BVI %%'
                AND phot.objid = frames.objid
                AND filter = '%(filtercol)s'
                AND good = True
                AND circle(phot.coord,0) <@ circle(point(%(ra).11f,%(dec).11f), 0.6/3600.0)
                AND flags<5;""" % params
            result = self.wifsip.query(query)
            #print result
            if len(result) == 0:
                continue
            elif len(result)==1:
                starid, mags, corr, err = result[0]
                mag = mags - corr
                std = err
                n = 1
            elif len(result)>1:
                mags = np.array([r[1] for r in result])
                corrs = np.array([r[2] for r in result])
                errs = np.array([r[3] for r in result])
                starid = result[0][0]
                mags -= corrs
                std0 = np.std(mags)
                # remove outliers
                i = np.where(abs(mags-np.mean(mags))<std0)
                mags = mags[i]
                errs = errs[i]
                std = np.std(mags)
                
                try:
                    mag = np.average(mags, weights = 1./errs)
                except TypeError:
                    print mags, corrs, errs
                    continue
                except ZeroDivisionError:
                    print mags, corrs, errs
                    mag = np.mean(mags)
                std = np.std(mags)
                n = len(mags)
            params['starid'] = starid
            params['mag'] = mag
            params['std'] = std
            params['n'] = n
            filterletter = filtercol.lower() 
            params['magfield'] = filterletter +'mag'
            params['errfield'] = filterletter +'mag_err'
            params['numfield'] = 'n' + filterletter
                
            print '%(starid)-25s, %(ra).6f, %(dec).6f, %(mag).4f, %(std).4f, %(n)3d' % params
            query = """INSERT INTO ngc6633 
            (starid,ra,dec,%(magfield)s,%(errfield)s,%(numfield)s)
            VALUES ('%(starid)s', %(ra).11f, %(dec).11f, %(mag)f, %(std)f, %(n)d);
            """ % params
            query = query.replace('nan', 'NULL')
            try:
                self.wifsip.execute(query)
            except IntegrityError:
                continue
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='NGC 6633 catalogue builder')
    parser.add_argument('-filter', help='filter color', default='V', dest='filtercol')
    parser.add_argument('inputfile', help='input file')

    args = parser.parse_args()
                
    nc = NGC6633cat()
    nc.load_positions(filename = args.inputfile)
    nc.build(filtercol = args.filtercol)