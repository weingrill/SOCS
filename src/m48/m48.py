#!/usr/bin/python
'''
Created on May 8, 2013

@author: jwe <jweingrill@aip.de>

Data reduction Class for M48 observation
'''

import logging
from datasource import DataSource

def autocorrelate(t, y):
    """
    calculate the autocorrelation on irregular data
    """
    import numpy as np
    from scipy.interpolate import interp1d

    #determine the median time step dt:
    dt = np.median(t-np.roll(t, 1))
    t0 = min(t)
    t1 = max(t)
    n = int((t1-t0)/dt)
    nt = np.linspace(t0, t1, n)
    f = interp1d(t, y, kind='linear')
    fn = f(nt)
    
    m = np.mean(fn)
    # we do padding
    k = 2
    x = fn - m
    # complex FFT, since we need the conjugate
    s = np.fft.fft(x, k*n)
    ac = np.real(np.fft.ifft(s*s.conjugate()))
    ac /= ac[0]
    lag = nt-nt[0]
    return ac[:n], lag[:n]

class Ngc1647(object):
    '''
    classdocs
    '''


    def __init__(self, path):
        '''
        Constructor
        '''
        self.stars = []
        self.path = path
        self.age = 10**8.158/1e6 # in Gyr from Webda
        self.ebv = 0.37 # from Webda
        self.dm = 9.81 # from Webda
        logging.basicConfig(filename='m48_analysis.log', 
                            format='%(asctime)s %(message)s')
    
        
    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        query = "select id from ngc1647stars where period is NULL;"
        #query = """SELECT id 
        #FROM ngc1647stars 
        #WHERE vmag<4.762*bv + 10.4762
        #AND period is NULL;"""
        logging.info('fetching stars ...')
        result = wifsip.query(query)
        logging.info('... %d stars found' % len(result))
        self.stars = [s[0] for s in result]
    
    def lightcurve_fromdb(self, starid):
        """
        extract a single lightcurve from the database
        and return epoch (hjd), magnitude and error
        """ 
        from datasource import DataSource
        import numpy as np
        
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT frames.hjd, phot.mag_auto, phot.magerr_auto, phot.flags
                FROM frames, matched, phot
                WHERE id LIKE '%s'
                AND filter LIKE 'rp'
                AND frames.objid = matched.objid
                AND (phot.objid,phot.star) = (matched.objid,matched.star)
                ORDER BY hjd;""" % (starid)
    # AND hjd>2455473.5 AND hjd<2455477 AND frames.good
        logging.info('fetching star %s' % starid)        
    # AND frames.good
    #          AND hjd>2455470 AND hjd<2455510
        data = wifsip.query(query)
        wifsip.close()
        hjd = np.array([d[0] for d in data])
        mag = np.array([d[1] for d in data])
        err = np.array([d[2] for d in data])
        return (hjd, mag, err)
    
    
    def periods(self):
        """
        calculate periods from the lightcurves
        """
        import pylab as plt
        import numpy as np
        
        self.means = np.empty(len(self.stars))
        self.stds = np.empty(len(self.stars))
        #f = open('/work1/jwe/NGC1647/sigma.tsv','wt')
        for star in self.stars: 
            i = self.stars.index(star)
            t, m, e = self.lightcurve_fromdb(star)
            t -= min(t)
            mean = np.mean(m)
            std = np.std(m)
            self.means[i] = mean
            self.stds[i] = std
            print '%s: %.3f %.3f' % (star,mean, std)
            #f.write('%s\t%.3f\t%.3f\n' % (star,mean, std))
            plt.scatter(mean, std)
            ac, lag = autocorrelate(t,m)
            plt.subplot(211)
            plt.title(star)
            plt.hlines(mean,min(t),max(t),linestyle='--')
            plt.hlines(mean-std*3,min(t),max(t),linestyle='dotted')
            plt.hlines(mean+std*3,min(t),max(t),linestyle='dotted')
            plt.ylim(max(m),min(m))
            plt.scatter(t,m)
            plt.subplot(212)
            plt.plot(lag, ac)
            plt.vlines(1.,min(ac),max(ac), linestyle='dotted')
            plt.grid()
            plt.hlines(0,min(lag), max(lag),linestyle='--')
            plt.show()
        #f.close()
        plt.ylabel(r'$\sigma$')
        plt.xlabel('r\' mag ')
        plt.yscale('log')
        plt.show()

    def analysis(self):
        import numpy as np
        import matplotlib.pylab as mpl
        from PyAstronomy.pyTiming import pyPDM

        for star in self.stars: 
            #i = self.stars.index(star)
            print star,'\t',
            t, m, e = self.lightcurve_fromdb(star)
            t -= min(t)
                
        
            S = pyPDM.Scanner(minVal=0.1, maxVal=10.0, dVal=10./86400., mode="period")
            P = pyPDM.PyPDM(t, m)
            p1, t1 = P.pdmEquiBinCover(100, 3, S)
            # For comparison, carry out PDM analysis using 10 bins equidistant
            # bins (no covers).
            #p2, t2 = P.pdmEquiBin(100, S)
        
            print '%.3f\t%.2f' % (p1[np.argmin(t1)], min(t1)) 
            # Show the result
            mpl.figure(facecolor='white')
            mpl.title("Result of PDM analysis")
            mpl.xlabel("Period")
            mpl.ylabel("Theta")
            mpl.plot(p1, t1)
            #mpl.plot(p2, t2, 'gp-')
            #mpl.legend(["pdmEquiBinCover", "pdmEquiBin"])
            #mpl.savefig('/work1/jwe/NGC1647/results/'+star+'.png')
            mpl.show()       

    def store_pdm(self, star, periods, thetas):
        try:
            f = open('/work2/jwe/NGC1647/results/'+star+'.tsv', 'wt')
            for s in zip(periods,thetas):
                f.write('%f\t%f\n' % s)
        finally:
            f.close()
    
    def store(self, star, period=None, theta=None):
        from datasource import DataSource
     
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        query = "UPDATE ngc1647stars"
        if not period is None:
            query += " SET period=%f" % period
        if not theta is None:
            query += ",  theta=%f" % theta
        query += " WHERE id like '%s';" %star
        try:
            wifsip.execute(query)
        finally:
            wifsip.close()
        
            
    def my_analysis(self):
        import numpy as np
        from pdm import pdm
        
        for star in self.stars: 
            #i = self.stars.index(star)
            t, m, _ = self.lightcurve_fromdb(star)
            try:
                t -= min(t)
            except ValueError:
                comment = '%s\t no data' % star
            # look at ten days or at most at the length of dataset
            else:
                length = min([max(t), 10.0])
                try:
                    p1, t1 = pdm(t, m, 0.1, length, 1./86400.)
                    period = p1[np.argmin(t1)]
                    theta = min(t1)
                    comment = '%s\t%.3f\t%.2f' % (star, period, theta)
        
                    self.store(star, period=period, theta=theta) 
                    self.store_pdm(star, p1, t1)
                except (ValueError,ZeroDivisionError):
                    comment = '%s\t no period' % star
            logging.info( comment)
            
    def make_cmd(self):
        import pylab as plt
        import numpy as np
        from datasource import DataSource
    
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = "SELECT vmag, bv FROM ngc1647stars;"
        
        data = wifsip.query(query)
        wifsip.close()
        vmag = np.array([d[0] for d in data])
        bv = np.array([d[1] for d in data])
        plt.scatter(bv,vmag, edgecolor='none', alpha=0.75)
        
        k = 8/(2.0-0.32)
        d = 20 - 2.0*k
        x = np.linspace(0.0, 2.5, 10)
        y = k*x+d
        plt.plot(x, y, linestyle='dashed', color='k')
        plt.ylim(21.0, 8.0)
        plt.xlim(0.0,2.2)
        plt.xlabel('B - V')
        plt.ylabel('V [mag]')
        plt.grid()
        plt.savefig('/work1/jwe/Dropbox/NGC1647/plots/ngc1647cmd.pdf')
        #plt.show()
        plt.close()

    def make_cpd(self):
        import pylab as plt
        import numpy as np
        from datasource import DataSource
    
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT bv, period, theta 
                    FROM ngc1647stars 
                    WHERE vmag>4.762*bv + 10.4762 and theta>0.5;"""
        data = wifsip.query(query)
        
        bv = np.array([d[0] for d in data])
        period = np.array([d[1] for d in data])
        theta = np.array([d[2] for d in data])

        query = """SELECT bv, period, theta 
                    FROM ngc1647stars 
                    WHERE vmag<4.762*bv + 10.4762 
                    AND theta>0.5;"""
        data = wifsip.query(query)
        wifsip.close()
        bv_ms = np.array([d[0] for d in data])
        period_ms = np.array([d[1] for d in data])
        theta_ms = np.array([d[2] for d in data])

        
        import gyroage
        
        bv170 = np.linspace(0.5, 1.6, num=20)
        P = gyroage.gyroperiod(bv170, 170.0)
        P01 = gyroage.gyroperiod(bv170, 170.0, P0=0.1)
        P33 = gyroage.gyroperiod(bv170, 170.0, P0=3.3)
        
        plt.plot(bv170, P, color='r')
        plt.plot(bv170, P01, color='r', linestyle='dashed')
        plt.plot(bv170, P33, color='r', linestyle='dashed')
        
        
        plt.scatter(bv-self.ebv, period, s=(1.0-theta)*40., edgecolor='none',alpha=0.5)
        plt.scatter(bv_ms-self.ebv, period_ms, s=(1.0-theta_ms)*40., edgecolor='none',
                    alpha=0.75, facecolor='green')
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('period [days]')
        plt.ylim(0.0, 10.0)
        plt.xlim(0.0, 2.2)
        plt.grid()
        plt.savefig('/work1/jwe/Dropbox/NGC1647/plots/ngc1647cpd.pdf')
        plt.show()
        plt.close()
        
        
            
if __name__ == '__main__':
    import matplotlib
    
    ngc1647 =  Ngc1647('/work1/jwe/NGC1647/data/')
    ngc1647.getstars()
    #print ngc1647.stars
    #ngc1647.export_lightcurves()
    #ngc1647.periods()
    #ngc1647.analysis()
    #ngc1647.my_analysis()
    ngc1647.make_cpd()
