#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 5, 2014

@author: jwe
'''
try:
    import config
except ImportError:
    print "Please define config file first"
    exit()

class Photometry(object):
    '''
    Photometry class has the following tasks:
    
    * create the db table
    * optionally clear the values for recalculation
    * build up a list of stars
    * collect B and V magnitudes for each star
    '''


    def __init__(self, objname='', filtercol='V', dbname=''):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        
        if objname<>'': # 'M48 BVI'
            self.objname=objname
        else:
            raise(ValueError,'objname not set')
            
        
        if filtercol in ('B','V','R', 'I'):
            self.filtercol = filtercol
        else:
            raise(ValueError, 'unknown filter color')
        
        self.frames=[]
        if dbname<>'': # 'M48stars'
            self.dbname=dbname
        else:
            raise(ValueError,'tablename not set')
        
        print 'filter %s' % self.filtercol

    def createtable(self):
        '''
        create table for the photometry
        '''
        if not raw_input('press Y to erase '+self.dbname)=='Y':
            return
        
        query = """
        DROP TABLE IF EXISTS %(dbname)s;
        CREATE TABLE %(dbname)s(
         starid varchar(25),
         bv real,
         vmag real default 0,
         vmag_err real,
         bmag real default 0,
         bmag_err real,
         period real,
         period_err real,
         amp real,
         amp_err real,
         nv integer default 0,
         nb integer default 0,
         ra double precision,
         dec double precision,
         coord point,
         PRIMARY KEY (starid));
        GRANT SELECT ON %(dbname)s TO public;
        CREATE INDEX idx_%(dbname)s_coords ON %(dbname)s USING GIST (circle(coord,1.0/3600.));
        """ % {'dbname':self.dbname}
        self.wifsip.execute(query)
        print "table '%s' created" % self.dbname

    def cleartable(self):
        """
        clears photometric values from the table for recalculation
        """
        if not raw_input('press Y to clear '+self.dbname)=='Y':
            return
        query = """
        UPDATE %s
        SET vmag = 0, vmag_err = NULL, bmag = 0, bmag_err =NULL, nv = 0, nb = 0, bv = NULL;
        """ % self.dbname
        self.wifsip.execute(query)
        print "table '%s' cleared" % self.dbname
        
    def getframes(self):
        for field in ['C','NW','NE','SW','SE']:
            objname = self.objname+' '+field
            query = """SELECT object, objid, abs(corr)
             FROM frames
             WHERE object LIKE '%s'
             AND filter LIKE '%s'
             AND NOT corr IS NULL
             ORDER BY abs(corr)
             limit 5;""" % (objname, self.filtercol)
            
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

    def addstar(self, starid, mag, ra, dec):
        # identify by coordinates, if the star is already in table
        query = """SELECT starid
        FROM %(dbname)s
        WHERE circle(point(%(ra)f,%(dec)f),0)<@circle(coord,1.0/3600.)
        ORDER BY point(%(ra)f,%(dec)f)<->coord
        LIMIT 1;""" % {'dbname': self.dbname, 'ra':ra, 'dec':dec}
        result = self.wifsip.query(query)
        oldstar = 0
        # if not: append new star
        if self.filtercol=='B': mname, nname= ('bmag','nb')
        elif self.filtercol=='V': mname, nname=( 'vmag','nv')

        if len(result)==0:
            oldstar = 0
            query="""INSERT INTO %s (starid, %s, %s, ra, dec, coord)
            VALUES ('%s', %f, 1, %f, %f, point(%f,%f))""" % \
            (self.dbname, mname, nname, starid, mag, ra, dec, ra, dec)
             
        # if star exists: add up magnitudes, increase counter
        else:
            oldstar = 1
            oldid = result[0][0]
            query = """UPDATE %(dbname)s
            SET %(mname)s = %(mname)s + %(mag)f, %(nname)s = %(nname)s + 1
            WHERE starid = '%(oldid)s';
            """ % {'dbname': self.dbname, 'mag':mag, 'oldid':oldid}
        self.wifsip.execute(query)
        return oldstar

    def update_magnitudes(self):
        if self.filtercol=='B': magfield, nfield = 'bmag','nb'
        elif self.filtercol=='V':magfield, nfield = 'vmag','nv'
        query = """UPDATE %(dbname)s
        SET %(magfield)s=%(magfield)s/%(nfield)s
        WHERE %(nfield)s>1;
        UPDATE %(dbname)s
        SET %(magfield)s=NULL
        WHERE %(nfield)s=0;""" % \
        {'dbname':self.dbname, 'magfield':magfield, 'nfield':nfield}
         
        self.wifsip.execute(query)   
    
    def update_sigmas(self):
        """
        update the photometric erros after the magnitude has been calculated
        """
        import numpy as np
        
        if self.filtercol=='V': field= 'vmag'
        elif self.filtercol=='B': field= 'bmag'
        query = """SELECT starid, coord 
            FROM %s 
            WHERE (NOT bv IS NULL) AND (%s_err IS NULL);""" % (self.dbname,field)
        starlist = self.wifsip.query(query)
        for star in starlist:
            print '%5d '% starlist.index(star),
            print '%-24s: %-25s' % star,
            query = """SELECT phot.objid, mag_auto-corr 
                FROM phot, frames
                WHERE object like '%s %%'
                AND phot.objid=frames.objid
                AND filter='%s'
                AND flags<8
                AND point%s <@ circle(phot.coord,1./3600.)
                ORDER BY abs(corr)
                LIMIT 5;""" % (self.objname, self.filtercol,star[1])
            result = self.wifsip.query(query)
            mags= np.array([r[1] for r in result ])
            try:
                err = np.std(mags)
                print mags,
                print '%.3f %.4f' % (np.mean(mags),err) 
                if np.isfinite(err):
                    query = "UPDATE %s SET %s_err=%f WHERE starid='%s';" % \
                    (self.dbname, field, err, star[0]) 
                    self.wifsip.execute(query)
            except TypeError:
                print 'no data'
    
    def update_bv(self):
        """
        just calculate the B-V
        """
        query = "UPDATE %s SET bv = bmag-vmag;" % self.dbname
        self.wifsip.execute(query)

def make_cmd(dbname):
    import pylab as plt
    import numpy as np
    from datasource import DataSource

    wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
    query = """UPDATE %s 
               SET bv = (bmag-vmag)+(0.761837-0.060177);""" % dbname
    
    wifsip.execute(query)
    
    query = """SELECT vmag+0.060177, bv 
                FROM %s 
                WHERE NOT bv is NULL and nv>1 and nb>1;""" % dbname
    
    data = wifsip.query(query)
    wifsip.close()
    vmag = np.array([d[0] for d in data])
    bv = np.array([d[1] for d in data])
    plt.scatter(bv,vmag, edgecolor='none', alpha=0.75, s=4, c='k')
    
    plt.ylim(22.0, 7.0)
    plt.xlim(-0.2, 2.0)
    plt.xlabel('B - V')
    plt.ylabel('V [mag]')
    plt.title(dbname)
    plt.grid()
    plt.savefig('/work2/jwe/m48/m48cmd.pdf')
    #plt.show()
    plt.close()


phot = Photometry()
phot.createtable()
phot.cleartable()
phot.getframes()
phot.filltable('2005')

