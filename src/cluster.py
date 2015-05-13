#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Cluster(dict):
    '''
    Cluster object to access the clusters table on wifsip database
    '''
    def __init__(self, name):
        '''
        Constructor:
        set the name of the cluster and open database
        '''
        from datasource import DataSource
        
        self.name = name
        
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        
    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'clusters';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        query = """SELECT * from clusters where name = '%s'""" % self.name
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        return values

    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM clusters 
        WHERE name='%s';""" % (key, self.name))
        return result[0][0]  
    
    @property
    def coordinates(self):
        """
        return the coordinates as a tuple of floats
        """
        return (float(self['ra']), float(self['dec']))  

    @property
    def coordinatestring(self):
        """
        return the coordinates as strings
        we don't need a high precision for the cluster coordinates
        """
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        
        ra = self['ra']
        dec = self['dec']
        c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)  # @UndefinedVariable
        ra_str =  c.ra.to_string(unit=u.hourangle,sep=' ', precision=1)  # @UndefinedVariable
        dec_str =  c.dec.to_string(sep=' ', precision=0)
        return ra_str, dec_str
    
    @property
    def distancemodulus(self):
        """
        calculates the distance modulus of the cluster from the distance in 
        parsec
        """
        from numpy import log10
        return 5.0*log10(self['d']) - 5.0
        


if __name__ == '__main__':
    c = Cluster('NGC 6709')
    print c['ra'],c['dec']
    print c['d'],c['diam']
    print c.coordinates
    print c.coordinatestring 
    print c.distancemodulus       