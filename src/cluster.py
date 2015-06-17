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
        self.separator = ' '
        self.ra_precision = 1
        self.dec_precision = 0
        self._keys = self.keys()
        self._values = self.values()
        for k, v in zip(self._keys, self._values):
            self[k] = v
        
        
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
        ra_str =  c.ra.to_string(unit=u.hourangle, # @UndefinedVariable
                                 sep=self.separator, 
                                 precision=self.ra_precision)  
        dec_str =  c.dec.to_string(sep=self.separator, 
                                   precision=self.dec_precision)
        return ra_str, dec_str
    
    @property
    def distancemodulus(self):
        """
        calculates the distance modulus of the cluster from the distance in 
        parsec
        """
        from numpy import log10
        return 5.0*log10(self['d']) - 5.0
        
    @property
    def solarmagnitude(self):
        """
        calculates the magnitude of a solar like star at the cluster distance
        """
        from numpy import log10
        return 5.0*log10(self['d']) - 5.0 + 4.862

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Cluster table query')
    parser.add_argument('clustername', help='name of the cluster')
    parser.add_argument('-k', '--keys', action='store_true',
                        help='list of available keys')
    parser.add_argument('-p', help='return specific parameter')
    
    args = parser.parse_args()
    
    c = Cluster(args.clustername)
    c.separator = ':'
    if args.keys:
        for key in c.keys(): print key 
    elif args.p:
        print c[args.p]
    else:
        print 'cluster:          %s' % c['name']
        print 'coordinates:      %s %s' % c.coordinatestring
        print 'proper motions    %.2f %.2f' % (c['pmra'],c['pmdec'])
        print 'radial velocity   %.2f' % c['rv']
        print 'diameter:         %.1f' % c['diam']
        print 'distance:         %d' % c['d']
        print 'E(B - V):         %.2f' % c['ebv']
        try:
            print '[Fe/H]:           %.2f' % c['me']
        except TypeError:
            pass
        print 'distance modulus: %.2f' % c.distancemodulus 
        print 'solar magnitude:  %.2f' % c.solarmagnitude
        print 'log age           %.3f' % c['logage']
            