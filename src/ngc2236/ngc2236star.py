'''
Created on Jul 18, 2014

@author: jwe
'''

import config
import logging
logging.basicConfig(filename='/work2/jwe/NGC2236/ngc2236star.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC2236 analysis')

class LightCurve(object):
    """
    Lightcurve object for NGC2236
    fetch lightcurve either from db or from file.
    ability to store lightcurve to file
    """
    def __init__(self, corotid):
        from datasource import DataSource
        self.corotid = corotid
        self.corot = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.fromdb()


    def fromfile(self, filename=None):
        """
        loads the lightcurve from a file.
        if the filename is not given it is assembled from the starid
        """
        import numpy as np
        
        if filename is None:
            filename = config.lightcurvespath+self.starid+'.dat'
        logger.info('load file %s' % filename)
        data = np.loadtxt(filename)
        
        self.hjd = data[:,0]
        self.mag = data[:,1]
        self.err = data[:,2]
                
        logger.info('%d datapoints' % len(self.hjd))
        
        return (self.hjd, self.mag, self.err)

class NGC2236Star(dict):
    '''
    class that interfaces the ngc2236 table on wifsip database
    '''
    def __init__(self, starid):
        from datasource import DataSource
        self.starid = starid
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        if '%' in self.starid:
            self.starid = self['starid']

    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'ngc2236';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        if '%' in self.starid:
            query = """SELECT * from ngc2236 where starid like '%s'""" % self.starid
        else:
            query = """SELECT * from ngc2236 where starid = '%s'""" % self.starid
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        if '%' in self.starid:
            self.starid=values[0]
        return values

    def __setitem__(self, key, value):
        if value is None:
            query = """UPDATE ngc2236 
            SET %s=NULL 
            WHERE starid='%s';""" % (key, self.starid)
        else:
            if type(value) is str:
                value = "'%s'" % value            
            query = """UPDATE ngc2236 
            SET %s=%s 
            WHERE starid='%s';""" % (key, str(value), self.starid)
        self.wifsip.execute(query)
        
    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM ngc2236 
        WHERE starid = '%s';""" % (key, self.starid))
        return result[0][0]    

    def lightcurve(self):
        lc = LightCurve(self.starid)
        return lc.hjd, lc.mag, lc.err