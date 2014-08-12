'''
Created on Jul 18, 2014

@author: jwe
'''

import logging
logging.basicConfig(filename='/work2/jwe/m48/m48_analysis.log', 
                    format='%(asctime)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')

class M48Star(dict):
    '''
    class that interfaces the m48stars table on wifsip database
    '''
    def __init__(self, starid):
        from datasource import DataSource
        self.starid = starid
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        if '%' in self.starid:
            self.starid = self['starid']

    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'm48stars';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        if '%' in self.starid:
            query = """SELECT * from m48stars where starid like '%s'""" % self.starid
        else:
            query = """SELECT * from m48stars where starid = '%s'""" % self.starid
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        if '%' in self.starid:
            self.starid=values[0]
        return values

    def __setitem__(self, key, value):
        if value is None:
            query = """UPDATE m48stars 
            SET %s=NULL 
            WHERE starid='%s';""" % (key, self.starid)
        else:            
            query = """UPDATE m48stars 
            SET %s=%s 
            WHERE starid='%s';""" % (key, str(value), self.starid)
        self.wifsip.execute(query)
        
    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM m48stars 
        WHERE starid like '%s';""" % (key, self.starid))
        return result[0][0]    
        
    def lightcurve(self):
        """
        extract a single lightcurve from the database
        and return epoch (hjd), magnitude and error
        """
        import numpy as np
        
        objid, star = self.starid.split('#')
        query = """SELECT id 
         FROM matched
         WHERE (matched.objid, matched.star) = ('%s','%s')
        """ % (objid, star)
        try:
            mid = self.wifsip.query(query)[0][0]
        except IndexError:
            logger.warning('no match found for starid %s' % (self.starid))
            return
        
        logger.info('fetching starid %s = %s' % (self.starid, mid))
        
        query = """SELECT frames.hjd, phot.mag_auto-corr, corr/2.0
                FROM frames, matched, phot
                WHERE matched.id LIKE '%s'
                AND frames.object like 'M 48 rot%%'
                AND filter LIKE 'V'
                AND NOT corr IS NULL
                AND frames.objid = matched.objid
                AND (phot.objid,phot.star) = (matched.objid,matched.star)
                AND phot.flags<4
                ORDER BY hjd;""" % (mid)
                
    
        data = self.wifsip.query(query)
        if len(data)<10:
            logger.error('insufficient data (%d) found for star %s' % (len(data),mid))
            return
        hjd = np.array([d[0] for d in data])
        mag = np.array([d[1] for d in data])
        err = np.array([d[2] for d in data])
                
        logger.info('%d datapoints' % len(hjd))
        
        return (hjd, mag, err)
