#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 22, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import config

class PhotDatabase(object):
    '''
    classdocs
    '''


    def __init__(self, filtercol):
        '''
        Constructor
        '''
        self.data = []
        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost, dictcursor=True)
        if not filtercol in ['u', 'v', 'b', 'y', 'hbn', 'hbw', 'I']:
            raise(ValueError)
        self.__filtercol=filtercol
        print self.__filtercol
        
        
    def fromfile(self, filename=None):
        import pyfits
        
        hdulist = pyfits.open(filename)
        self.data = hdulist[1].data
        self.keys = []
        for col in hdulist[1].columns:
            self.keys.append(col.name)
        hdulist.close()
    
    def todatabase(self):
        for d in self.data:
            record = {}
            for key in self.keys:
                record[key] = d[key]
            record['filter'] = self.__filtercol
            if self.__filtercol == 'I':
                record['filter'] = 'imag'
            query = """UPDATE ngc2236 
            SET %(filter)s = %(MAG_ISOCOR)f, %(filter)s_err = %(MAGERR_ISO)f 
            WHERE circle(coord,0) <@ circle(point(%(ALPHAWIN_J2000).11f,%(DELTAWIN_J2000).11f), 0.6/3600.0);""" % record
            self.wifsip.execute(query, commit=False)
            
        self.wifsip.commit()        

if __name__ == '__main__':
    for col in ['I']:#['u', 'v', 'b', 'y', 'hbn', 'hbw']:
        phdb = PhotDatabase(col)
        filename = '/work2/jwe/SOCS/NGC2236/frames/%s/NGC2236%s.fit' % (col,col)
        phdb.fromfile(filename)
        phdb.todatabase()