#!/usr/bin/python
'''
Created on May 8, 2013

@author: jwe <jweingrill@aip.de>

Data reduction Class for M48 observation
'''

import logging
#import matplotlib

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

class M48(object):
    '''
    classdocs
    '''


    def __init__(self, path):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.stars = []
        self.path = path
        self.age = 10**8.557/1e6 # in Gyr from Webda
        self.ebv = 0.031 # from Webda
        self.dm = 9.53 # from Webda
        logging.basicConfig(filename='m48_analysis.log', 
                            format='%(asctime)s %(message)s',
                            level=logging.INFO)
        
    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        
        query = """SELECT starid 
        FROM m48stars 
        WHERE period IS NULL 
        AND NOT bv IS NULL;"""
        
        logging.info('fetching stars ...')
        result = self.wifsip.query(query)
        logging.info('... %d stars found' % len(result))
        self.stars = [s[0] for s in result]
    
    def lightcurve_fromdb(self, starid):
        """
        extract a single lightcurve from the database
        and return epoch (hjd), magnitude and error
        """
        import numpy as np
        
        objid, star = starid.split('#')
        query = """SELECT id 
         FROM matched
         WHERE (matched.objid, matched.star) = ('%s','%s')
        """ % (objid, star)
        try:
            mid = self.wifsip.query(query)[0][0]
        except IndexError:
            logging.warning('no match found for starid %s' % (starid))
            self.stars.remove(starid)
            return
        
        logging.info('fetching starid %s = %s' % (starid, mid))
        
        query = """SELECT frames.hjd, phot.mag_auto, phot.magerr_auto
                FROM frames, matched, phot
                WHERE matched.id LIKE '%s'
                AND frames.object like 'M 48 rot%%'
                AND filter LIKE 'V'
                AND frames.objid = matched.objid
                AND (phot.objid,phot.star) = (matched.objid,matched.star)
                AND phot.flags<8
                ORDER BY hjd;""" % (mid)
                
    # AND frames.good
        data = self.wifsip.query(query)
        if len(data)<3:
            logging.error('no data found for star %s' % mid)
            return
        hjd = np.array([d[0] for d in data])
        mag = np.array([d[1] for d in data])
        err = np.array([d[2] for d in data])
        logging.info('%d datapoints' % len(hjd))
        return (hjd, mag, err)
    
    def periods(self):
        """
        calculate periods from the lightcurves
        """
        import pylab as plt
        import numpy as np
        
        self.means = np.empty(len(self.stars))
        self.stds = np.empty(len(self.stars))
        f = open('/work2/jwe/m48/sigma.tsv','wt')
        for starid in self.stars: 
            i = self.stars.index(starid)
            try:
                t, m, e = self.lightcurve_fromdb(starid)
            except TypeError:
                pass
            else:
                t -= min(t)
                mean = np.mean(m)
                std = np.std(m)
                self.means[i] = mean
                self.stds[i] = std
                print '%s: %.3f %.3f' % (starid, mean, std)
                f.write('%s\t%.3f\t%.3f\n' % (starid, mean, std))
                plt.scatter(mean, std)
                ac, lag = autocorrelate(t,m)
                plt.subplot(211)
                plt.title(starid)
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
                plt.savefig('/work2/jwe/m48/plots/%s.png' % starid)
                plt.close()
        f.close()
        

    def store_pdm(self, star, periods, thetas):
        try:
            f = open('/work2/jwe/m48/results/'+star+'.tsv', 'wt')
            for s in zip(periods,thetas):
                f.write('%f\t%f\n' % s)
        finally:
            f.close()
    
    def store(self, star, period=None, theta=None):
        query = "UPDATE m48stars"
        if not period is None:
            query += " SET period=%f" % period
        if not theta is None:
            query += ",  theta=%f" % theta
        query += " WHERE id like '%s';" %star
        self.wifsip.execute(query)
                    
    def analysis(self):
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
                length = min([max(t), 20.0])
                try:
                    p1, t1 = pdm(t, m, 0.1, length, 60./86400.)
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
        query = "SELECT vmag, bv FROM m48stars;"
        
        data = self.wifsip.query(query)
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
        plt.savefig('/work2/jwe/m48/plots/ngc1647cmd.pdf')
        #plt.show()
        plt.close()

    def make_cpd(self):
        import pylab as plt
        import numpy as np
        query = """SELECT bv, period, theta 
                    FROM m48stars 
                    WHERE vmag>4.762*bv + 10.4762 and theta>0.5;"""
        data = self.wifsip.query(query)
        
        bv = np.array([d[0] for d in data])
        period = np.array([d[1] for d in data])
        theta = np.array([d[2] for d in data])

        query = """SELECT bv, period, theta 
                    FROM m48stars 
                    WHERE vmag<4.762*bv + 10.4762 
                    AND theta>0.5;"""
        data = self.wifsip.query(query)
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
        plt.savefig('/work2/jwe/m48/plots/m48cpd.pdf')
        plt.show()
        plt.close()
        
    def __exit__(self):
        self.wifsip.close()
            
if __name__ == '__main__':
    
    m48 =  M48('/work2/jwe/m48/data/')
    logging.info('getting stars ...')
    m48.getstars()
    print '\n'.join(m48.stars)
    logging.info('periodsearch ...')
    m48.periods()
    m48.analysis()
    #ngc1647.make_cpd()
