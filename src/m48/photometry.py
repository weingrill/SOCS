#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 5, 2014

@author: jwe
'''


class Photometry(object):
    '''
    classdocs
    '''


    def __init__(self, filtercol='V'):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        if filtercol in ('B','V'):
            self.filtercol = filtercol
        else:
            raise(ValueError)
        
        self.frames=[]
        print 'filter %s' % self.filtercol

    def createtable(self):
        '''
        create table for the photometry
        '''
        if not raw_input('press Y to erase m48stars')=='Y':
            return
        
        query = """
        DROP TABLE IF EXISTS m48stars;
        CREATE TABLE m48stars(
         starid varchar(25),
         bv real,
         vmag real default 0,
         vmag_err real,
         bmag real default 0,
         bmag_err real,
         period real,
         period_err real,
         theta real,
         amp real,
         amp_err real,
         nv integer default 0,
         nb integer default 0,
         ra double precision,
         dec double precision,
         coord point,
         PRIMARY KEY (starid));
        GRANT SELECT ON m48stars TO public;
        CREATE INDEX idx_m48stars_coords ON m48stars USING GIST (circle(coord,1.0/3600.));
        """
        self.wifsip.execute(query)
        print "table 'm48stars' created"

    def cleartable(self):
        if not raw_input('press Y to clear m48stars')=='Y':
            return
        query = """
        UPDATE m48stars
        SET vmag = 0, vmag_err = NULL, bmag = 0, bmag_err =NULL, nv = 0, nb = 0, bv = NULL;
        """
        self.wifsip.execute(query)
        print "table 'm48stars' cleared"
        
    def getframes(self):
        for field in ['C','NW','NE','SW','SE']:
            query = """SELECT object, objid, abs(corr)
             FROM frames
             WHERE object LIKE 'M 48 BVI %s'
             AND filter LIKE '%s'
             AND NOT corr IS NULL
             ORDER BY abs(corr)
             limit 5;""" % (field,self.filtercol)
            
            result = self.wifsip.query(query)
            if len(result)==0:
                print 'no frames found!'
            for r in result: 
                print '%s\t%s\t%.3f: ' % r
                objid = r[1]
                self.filltable(objid)
                self.frames.append(objid)
            self.wifsip.dropview('phot1')
        self.update_magnitudes()
        #print '\n'.join(self.frames)

    def filltable(self, objid):
        #get the stars from the phot table ...
        query = """
            SELECT phot.objid ||'#'|| star, mag_auto-corr, alphawin_j2000 , deltawin_j2000
            FROM phot,frames
            WHERE frames.objid='%s'
            AND phot.objid=frames.objid AND flags<8;""" % (objid)
        result = self.wifsip.query(query)
        #... and inject them into the m48stars 
        stars = len(result)
        if stars>400:
            print '%5d stars: ' % stars,
            oldstars = 0
            newstars = 0 
            for r in result:
                ostars = self.addstar(r[0], r[1], r[2], r[3])
                if ostars == 0:
                    newstars += 1
                else:
                    oldstars += 1 
            print '%5d old , %5d new' % (oldstars, newstars)
            self.wifsip.commit()
        else: print 'not enough stars (%d)' % stars
            #value = '%s\t%f\t%f\t%f' % r
            #print value
        
            

    def addstar(self, starid, mag, ra, dec):
        # identify by coordinates, if the star is already in table
        query = """SELECT starid
        FROM m48stars
        WHERE circle(point(%f,%f),0)<@circle(coord,1.0/3600.)
        ORDER BY point(%f,%f)<->coord
        LIMIT 1;""" % (ra, dec, ra, dec)
        result = self.wifsip.query(query)
        oldstar = 0
        # if not: append new star
        if len(result)==0:
            oldstar = 0
            if self.filtercol=='B':
                query="""INSERT INTO m48stars (starid, bmag, nb, ra, dec, coord)
                VALUES ('%s', %f, 1, %f, %f, point(%f,%f))""" % (starid, mag, ra, dec, ra, dec)
            elif self.filtercol=='V':
                query="""INSERT INTO m48stars (starid, vmag, nv, ra, dec, coord)
                VALUES ('%s', %f, 1, %f, %f, point(%f,%f))""" % (starid, mag, ra, dec, ra, dec)
            
        # if star exists: add up magnitudes, increase counter
        else:
            oldstar = 1
            oldid = result[0][0]
            if self.filtercol=='B':
                query = """UPDATE m48stars
                SET bmag = bmag + %f, nb = nb + 1
                WHERE starid = '%s';
                """ % (mag, oldid)
            elif self.filtercol=='V':
                query = """UPDATE m48stars
                SET vmag = vmag + %f, nv = nv + 1
                WHERE starid = '%s';
                """ % (mag, oldid)
        self.wifsip.execute(query)
        return oldstar
    
    def update_magnitudes(self):
        if self.filtercol=='B':
            query = """UPDATE m48stars
            SET bmag=bmag/nb
            WHERE nb>1;
            UPDATE m48stars
            SET bmag=NULL
            WHERE nb=0;"""
        elif self.filtercol=='V':
            query = """    
            UPDATE m48stars
            SET vmag=vmag/nv
            WHERE nv>1;
            UPDATE m48stars
            SET vmag=NULL
            WHERE nv=0;""" 
        self.wifsip.execute(query)   
    
    def update_sigmas(self):
        import numpy as np
        
        if self.filtercol=='V': field= 'vmag'
        elif self.filtercol=='B': field= 'bmag'
        query = """SELECT starid, coord 
            FROM m48stars 
            WHERE (NOT bv IS NULL) AND (%s_err IS NULL);""" % (field)
        starlist = self.wifsip.query(query)
        for star in starlist:
            print '%5d '% starlist.index(star),
            print '%-24s: %-25s' % star,
            query = """SELECT phot.objid, mag_auto-corr 
                FROM phot, frames
                WHERE object like 'M 48 BVI %%'
                AND phot.objid=frames.objid
                AND filter='%s'
                AND flags<8
                AND point%s <@ circle(phot.coord,1./3600.)
                ORDER BY abs(corr)
                LIMIT 5;""" % (self.filtercol,star[1])
            result = self.wifsip.query(query)
            mags= np.array([r[1] for r in result ])
            try:
                err = np.std(mags)
                print mags,
                print '%.3f %.4f' % (np.mean(mags),err) 
                if np.isfinite(err):
                    query = "UPDATE m48stars SET %s_err=%f WHERE starid='%s';" % (field, err, star[0]) 
                    self.wifsip.execute(query)
            except TypeError:
                print 'no data'
    
    def update_bv(self):
        query = "UPDATE m48stars SET bv = bmag-vmag;"
        self.wifsip.execute(query)
        
def make_cmd():
    import pylab as plt
    import numpy as np
    from datasource import DataSource

    wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
    query = "UPDATE m48stars SET bv = (bmag-vmag)+(0.761837-0.060177);"
    
    wifsip.execute(query)
    query = "SELECT vmag+0.060177, bv FROM m48stars WHERE NOT bv is NULL and nv>1 and nb>1;"
    
    data = wifsip.query(query)
    wifsip.close()
    vmag = np.array([d[0] for d in data])
    bv = np.array([d[1] for d in data])
    plt.scatter(bv,vmag, edgecolor='none', alpha=0.75, s=4, c='k')
    
    plt.ylim(22.0, 7.0)
    plt.xlim(-0.2, 2.0)
    plt.xlabel('B - V')
    plt.ylabel('V [mag]')
    plt.title('M 48')
    plt.grid()
    plt.savefig('/work2/jwe/m48/m48cmd.pdf')
    #plt.show()
    plt.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='M48 photometry')
    parser.add_argument('--create', action='store_true', help='create table')
    parser.add_argument('--clear', action='store_true', help='clear table')
    parser.add_argument('--getframes', action='store_true', help='get frames')
    parser.add_argument('--sigmas', action='store_true', help='plot cmd')
    parser.add_argument('--bv', action='store_true', help='update B-V')
    parser.add_argument('filter', metavar='V', type=str, help='filter color to process')
    
    args = parser.parse_args()
    
    phot = Photometry(filtercol=args.filter)
    if args.create: phot.createtable()
    if args.clear: phot.cleartable()
    if args.getframes: phot.getframes()
    if args.sigmas: phot.update_sigmas()
    if args.sigmas: phot.update_bv()
