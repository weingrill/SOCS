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

    def createtable(self):
        '''
        create table for the photometry
        '''
        
        query = """
        DROP TABLE IF EXISTS m48stars;
        CREATE TABLE m48stars(
         starid varchar(25),
         bv real,
         vmag real,
         vmag_err real,
         bmag real,
         bmag_err real,
         period real,
         period_err real,
         theta real,
         amp real,
         amp_err real,
         nv integer,
         nb integer,
         ra double precision,
         dec double precision,
         coord point,
         PRIMARY KEY (starid));
        GRANT SELECT ON m48stars TO public;
        CREATE INDEX idx_m48stars_coords ON m48stars USING GIST (circle(coord,3.0/3600.));
        """
        self.wifsip.execute(query)
        print "table 'm48stars' created"

    def getframes(self):
        for field in ['C']: #['C','NW','NE','SW','SE']:
            query = """SELECT object, objid, abs(corr)
             FROM frames
             WHERE object LIKE 'M 48 BVI %s'
             AND filter LIKE '%s'
             AND NOT corr IS NULL
             AND good=True
             ORDER BY abs(corr)
             limit 5;""" % (field,self.filtercol)
            
            result = self.wifsip.query(query)
            if len(result)==0:
                print 'no frames found!'
            for r in result: 
                print '%s\t%s\t%.3f: ' % r,
                objid = r[1]
                self.filltable(objid)
                
                self.frames.append(objid)
            self.wifsip.dropview('phot1')
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
        if stars>800:
            print '%d stars' % stars
            for r in result:
                self.addstar(r[0], r[1], r[2], r[3])
        else: print 'not enough stars (%d)' % stars
            #value = '%s\t%f\t%f\t%f' % r
            #print value
            

    def addstar(self, starid, mag, ra, dec):
        # identify by coordinates, if the star is already in table
        query = """SELECT starid, (point(%f,%f)<->coord)*3600.
        FROM m48stars
        WHERE circle(point(%f,%f),0)<@circle(coord,3.0/3600.)
        ORDER BY point(%f,%f)<->coord
        LIMIT 1;""" % (ra, dec, ra, dec, ra, dec)
        result = self.wifsip.query(query)
        print result
        # if not: append new star
        if len(result)==0:
            print 'new'
            if self.filtercol=='B':
                query="""INSERT INTO m48stars (starid, bmag, nb, ra, dec, coord)
                VALUES ('%s', %f, 1, %f, %f, point(%f,%f))""" % (starid, mag, ra, dec, ra, dec)
            elif self.filtercol=='V':
                query="""INSERT INTO m48stars (starid, vmag, nv, ra, dec, coord)
                VALUES ('%s', %f, 1, %f, %f, point(%f,%f))""" % (starid, mag, ra, dec, ra, dec)
            self.wifsip.execute(query)
        # if star exists: add up magnitudes, increase counter
        else:
            oldid = result[0][0]
            print '= %s' % oldid
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
    
    def update_magnitudes(self):
        query = """UPDATE m48stars
        SET bmag=bmag/nb
        WHERE nb>1;
        UPDATE m48stars
        SET vmag=vmag/nv
        WHERE nv>1;
        UPDATE m48stars
        SET bv=bmag-vmag;""" 
        self.wifsip.execute(query)   

if __name__ == '__main__':
    p = Photometry(filtercol='V')
    p.createtable()
    p.getframes()
        
    p = Photometry(filtercol='B')
    p.getframes()
    p.update_magnitudes()