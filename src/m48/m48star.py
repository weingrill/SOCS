'''
Created on Jul 18, 2014

@author: jwe
'''

import config
import logging
logging.basicConfig(filename=config.projectpath+'m48star.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')

class LightCurve(object):
    """
    Lightcurve object for M48
    fetch lightcurve either from db or from file.
    ability to store lightcurve to file
    """
    def __init__(self, starid):
        self.starid = starid
        self.fromfile()

    def fromfile(self, filename=None):
        """
        loads the lightcurve from a file.
        if the filename is not given it is assembled from the starid
        """
        import numpy as np
        
        if filename is None:
            filename = config.lightcurvespath+self.starid+'.dat'
        logger.info('load file %s' % filename)
        self.hjd, self.mag, self.err = np.loadtxt(filename, unpack = True)
        
        logger.info('%d datapoints' % len(self.hjd))
        
        return (self.hjd, self.mag, self.err)

    def tofile(self, filename=None):
        """
        stores the lightcurve to a file
        """
        import numpy as np
        if filename is None:
            filename = config.lightcurvespath+self.starid+'.dat'
            
        a = np.column_stack((self.hjd,self.mag,self.err))
        np.savetxt(filename, a, fmt='%.6f %.3f %.4f')        

    def rebin(self, interval = 512./86400., medianbins=False):
        """
        rebin to new interval using the mean of each bin
        interval determines the length of each bin
        medianbins calculates the median in each bin otherwise the mean is taken
        """
        from numpy import seterr, zeros, isnan, compress, arange, mean
        data = self.hjd
        # ...+interval so that the last bin includes the last epoch
        bins = arange(self.hjd[0], self.hjd[-1]+interval, interval)
        nbins = len(bins)-1
        t = zeros(nbins)
        f = zeros(nbins)
        # adopted from Ian's Astro-Python Code v0.3
        # http://www.mpia-hd.mpg.de/homes/ianc/python/_modules/tools.html
        # def errxy()
        idx = [[data.searchsorted(bins[i]), \
                data.searchsorted(bins[i+1])] for i in range(nbins)]
        seterr(invalid='ignore')
        for i in range(nbins):
            f[i] = mean(self.mag[idx[i][0]:idx[i][1]])
            t[i] = mean(self.hjd[idx[i][0]:idx[i][1]])
                
        seterr(invalid='warn')
        valid = ~isnan(t)
        self.mag = compress(valid,f)
        self.hjd = compress(valid,t)
        
    def sigma_clip(self):
        from functions import sigma_clip
        self.hjd, self.mag = sigma_clip(self.hjd, self.mag)

    def detrend(self, degree = 1):
        """
        detrend the lightcurve with a polynomial of the order degree
        """
        from numpy import polyfit, polyval
        par = polyfit(self.hjd, self.mag, degree)
        self.mag -= polyval(par, self.hjd)

    @property
    def data(self):
        return self.hjd, self.mag
    

class M48Star(dict):
    '''
    class that interfaces the m48stars table on wifsip database
    '''
    def __init__(self, starid, tab= None):
        from datasource import DataSource
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        if starid is None:
            self.starid = self._staridfromtab(tab)
        else:
            self.starid = starid
        
        if '%' in self.starid:
            self.starid = self['starid']
    
    def _staridfromtab(self, tab):       
        result = self.wifsip.query("""SELECT starid 
        FROM m48stars 
        WHERE tab = %d;""" % tab)
        return result[0][0]

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
            if type(value) is str:
                value = "'%s'" % value            
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
        return LightCurve(self.starid)
    
    def cleanspectrum(self):
        #20140305A-0003-0017#359 <-- starid
        #       6A-0097-0013#1648.ncfile <-- filename
        #01234567
        from numpy import loadtxt
        filename = config.datapath+'clean/%s.ncfile' % self.starid[7:]
        a = loadtxt(filename)
        freq = a[:1600,0]
        amp = a[:1600,1]
        return freq, amp