'''
Created on Jul 18, 2014

@author: jwe
'''

import config
import logging
logging.basicConfig(filename=config.logpath+'ngc2236star.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC2236 analysis')

from lightcurve import LightCurve

class NGCLightCurve(LightCurve):
    """
    Lightcurve object for NGC2236
    fetch lightcurve either from db or from file.
    ability to store lightcurve to file
    """
    def _filename(self, corotid):
        """
        translates corotid to filename
        """
        from datasource import DataSource
        self.corotid = corotid
        self.corot = DataSource(database='corot', user='sro', host='pina.aip.de')
        
        query = """SELECT run_code, hlfccdid, win_id 
        FROM corot 
        WHERE corotid = %d;""" % self.corotid
        result = self.corot.query(query)
        
        par = {'run': result[0][0],
               'half': result[0][1].rstrip('RL'), 
               'win': result[0][2]}
        filename = '/work2/jwe/CoRoT/%(run)s/data/%(run)s_%(half)s_%(win)04d.fits' % par
        logger.info('%d = %s' % (corotid,filename))
        return filename
    
    def __init__(self, corotid, channel=None):
        """
        Constructor
        """
        from numpy import diff, where, nan
        import pyfits

        hdulist = pyfits.open(self._filename(corotid))
        hdulist.verify('silentfix')
        self.hdr = hdulist[0].header
        tbdata = hdulist[1].data
        self.time = tbdata.field('DATEHEL')
        status = tbdata.field('STATUS')
        if channel is None:
            if 'WHITEFLUX_IMAG' in tbdata.names:
                self.flux = tbdata.field('WHITEFLUX_IMAG')
            elif 'WHITEFLUX' in tbdata.names:
                self.flux = tbdata.field('WHITEFLUX')
        if channel in ['r','red']:
            if 'REDFLUX' in tbdata.names:
                self.flux = tbdata.field('REDFLUX')
        if channel in ['g','green']:
            if 'GREENFLUX' in tbdata.names:
                self.flux = tbdata.field('GREENFLUX')
        if channel in ['b','blue']:
            if 'BLUEFLUX' in tbdata.names:
                self.flux = tbdata.field('BLUEFLUX')
        
        hdulist.close()
        self.corotid = self.hdr['COROTID']

        i = where(status == 0)
        self.time = self.time[i] + 2451545.0 # 1.1.2000 12:00
        self.flux = self.flux[i]

        self.err  = abs(diff(self.flux)/diff(self.time))
        self.winid = self.hdr['RUN_CODE']+'_'+\
                     self.hdr['HLFCCDID'].rstrip('RL')+\
                     '_%04d' % self.hdr['WIN_ID']
        self.magr = 0.0
        
        if 'MAGNIT_R' in self.hdr:
            self.magr = self.hdr['MAGNIT_R']
        if 'MAGNIT_V' in self.hdr:
            self.magv = self.hdr['MAGNIT_V']
        else: self.magv = nan    
        
        if 'MAGNIT_B' in self.hdr:
            self.magb = self.hdr['MAGNIT_B']
        else: self.magb = nan
        if self.magb == 'NAN': self.magb = nan
        self.bv = float(self.magb)-float(self.magv)

class NGC2236Star(dict):
    '''
    class that interfaces the ngc2236 table on wifsip database
    '''
    def __init__(self, starid):
        """
        Constructor: connects to the database
        """
        from datasource import DataSource
        self.starid = starid
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        if '%' in self.starid:
            self.starid = self['starid']

    def keys(self):
        """
        returns the keys from the dict 
        """
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'ngc2236';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        """
        returns the values from the dict 
        """
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
        """
        setter method: sets value based on key in database table
        """
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
        """
        getter method: fetches value based on key from database table
        """
        result = self.wifsip.query("""SELECT %s 
        FROM ngc2236 
        WHERE starid = '%s';""" % (key, self.starid))
        return result[0][0]    

    def lightcurve(self):
        """
        return the lightcurve object
        """
        return NGCLightCurve(self['corotid'])
