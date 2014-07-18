'''
Created on Jul 18, 2014

@author: jwe
'''

import logging
logging.basicConfig(filename='/work2/jwe/m48/m48_analysis.log', 
                    format='%(asctime)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')

class M48Star(object):
    '''
    class that interfaces the m48stars table on wifsip database
    '''
    def __init__(self, starid):
        from datasource import DataSource
        self.starid = starid
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT bv, vmag, ra, dec, simbad 
        FROM m48stars
        WHERE starid='%s'""" % starid
        result = self.wifsip.query(query)[0]
        self.bv,self.vmag,self.ra,self.dec,self.simbad = result
    
    def _db_setvalue(self, param, value):
        if value is None:
            query = """UPDATE m48stars 
            SET %s=NULL 
            WHERE starid='%s';""" % (param, self.starid)
            self.wifsip.execute(query)
        else:            
            query = """UPDATE m48stars 
            SET %s=%f 
            WHERE starid='%s';""" % (param, value, self.starid)
        self.wifsip.execute(query)
                

    def _db_getvalue(self, param):
        result = self.wifsip.query("""SELECT %s 
        FROM m48stars 
        WHERE starid='%s';""" % (param, self.starid))
        return result[0][0]    
        
    @property
    def period(self):
        return self._db_getvalue('period')
        
    @period.setter
    def period(self, value):
        self._db_setvalue('period', value)

    @property
    def period_err(self):
        return self._db_getvalue('period_err')
        
    @period_err.setter
    def period_err(self, value):
        self._db_setvalue('period_err', value)


    @property
    def theta(self):
        return self._db_getvalue('theta')
        
    @theta.setter
    def theta(self, value):
        self._db_setvalue('theta', value)

    @property
    def amp(self):
        return self._db_getvalue('amp')

    @amp.setter
    def amp(self, value):
        #print value
        self._db_setvalue('amp', value)
     
    @property
    def amp_err(self):
        return self._db_getvalue('amp_err')

    @amp_err.setter
    def amp_err(self, value):
        self._db_setvalue('amp_err', value)

    @property
    def freq(self):
        return self._db_getvalue('freq')

    @freq.setter
    def freq(self, value):
        self._db_setvalue('freq', value)

    @property
    def s1(self):
        return self._db_getvalue('s1')

    @s1.setter
    def s1(self, value):
        self._db_setvalue('s1', value)

    @property
    def c1(self):
        return self._db_getvalue('c1')

    @c1.setter
    def c1(self, value):
        self._db_setvalue('c1', value)

    @property
    def s2(self):
        return self._db_getvalue('s2')

    @s2.setter
    def s2(self, value):
        self._db_setvalue('s2', value)

    @property
    def c2(self):
        return self._db_getvalue('c2')

    @c2.setter
    def c2(self, value):
        self._db_setvalue('c2', value)


        
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
        
        query = """SELECT frames.hjd, phot.mag_auto-corr, phot.magerr_auto
                FROM frames, matched, phot
                WHERE matched.id LIKE '%s'
                AND frames.object like 'M 48 rot%%'
                AND filter LIKE 'V'
                AND frames.good
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
