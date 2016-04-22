#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 21, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import config
import numpy as np

class NGC6633Table(object):
    def __init__(self):
        '''
        Constructor
        '''
        try:
            self._fromfile()
        except:
            self._fromdatabase()
            self._tofile()

    def _fromdatabase(self):
        """
        import the table from a database
        """
        from datasource import DataSource
        from astropy.coordinates import SkyCoord  # @UnresolvedImport
        from astropy import units as u
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        self.stars = []
        columns = self.wifsip.columns('ngc6633')
        
        params = {'columns': ', '.join(columns),
                  'key': columns[0]} 
        query = "SELECT %(columns)s FROM ngc6633 ORDER BY %(key)s;" % params 
        self._stars = self.wifsip.query(query)
        
        columns = self.wifsip.columns('ngc6633')
        data_types = self.wifsip.data_types('ngc6633')
        
        arraydata = []
        c = columns.index('coord')
        for star in self._stars:
            starlist = list(star)
            coord = star[c]
            ra, dec = coord.strip('()').split(',')
            
            newcoord = SkyCoord(ra=float(ra)*u.deg, dec=float(dec)*u.deg)
            starlist[c] = newcoord
            arraydata.append(tuple(starlist))
        
        dtype = zip(columns, data_types)
        #print dtype
        self.stars = np.array(arraydata, dtype = dtype)
        
    def _fromfile(self, filename=None):
        """
        unpickle the data from a file
        """
        import pickle
        
        if filename is None:
            filename = config.datapath+'/stars.pickle'
        with open(filename,'rb') as picklefile:
            self.stars = pickle.load(picklefile)
    
    def _tofile(self, filename=None):
        """
        pickle the data to a file
        """
        import pickle
        
        if filename is None:
            filename = config.datapath+'/stars.pickle'
        with open(filename, 'wb') as picklefile:
            pickle.dump(self.stars, picklefile)


if __name__ == '__main__':
    table = NGC6633Table()
    #print table.stars[:3]
    #print table.stars[0]['starid'],table.stars[0]['coord']
    #c = table.stars[0]['coord']
    #print type(c)
    
    import matplotlib.pyplot as plt
    
    plt.semilogy(table.stars[:]['vmag'], table.stars[:]['rms'],'ko', alpha=0.5)
    plt.show()

    
    #print a.stars[37]
