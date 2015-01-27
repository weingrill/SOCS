'''
Created on Jan 27, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from m48 import M48Analysis
import config
import logging
from m48star import M48Star            
import pylab as plt
import numpy as np
logging.basicConfig(filename=config.projectpath+'m48_analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('M48 analysis')

class M48ReAnalysis(M48Analysis):
    '''
    Inhertied class to redo the analysis for generating the lightcurve plots 
    without altering the table in the database.
    '''

    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        query = """SELECT starid 
        FROM m48stars 
        WHERE good
        ORDER BY tab;""" 

        logger.info('fetching stars ...')
        result = self.wifsip.query(query)
        logger.info('... %d stars found' % len(result))
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]


    def analysis(self, show=False):
        """perform a PDM analysis on each lightcurve"""
        from pdm import pdm
        from psd import ppsd
        from matplotlib import rcParams
        from functions import sigma_clip, phase
        print 'Analysis'

        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'font.size': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)

        minperiod = 1.2/24
        maxperiod = 15
        
        for starid in self.stars:
            star = M48Star(starid)
            print '%-24s '% starid,
            try:
                t, m, e = star.lightcurve()
                t -= min(t)
            except AttributeError:
                logger.error("Can't load lightcurve %s" % starid)
                print 'no lightcurve'
                continue
            
            if len(t)<50:
                logger.warn("%s: not enough datapoints" % starid)
                print 'not enough datapoints'
                continue                    
            # perform a 3sigma clipping
            self.t, self.m, self.e = sigma_clip(t, m, e)
            
            
            # perform a power spectrum analysis
            tpsa, mpsa = self.t, self.m- np.mean(self.m)
            n = len(tpsa)
            # zero padded lightcurves
            t_padded = np.zeros(4*n)
            t_padded[:n] = tpsa
            t_padded[n:] = np.linspace(max(tpsa),4*max(tpsa),3*n)
            m_padded = np.zeros(4*n)
            m_padded[:n] = mpsa
            
            px, f = ppsd(t_padded, 
                         m_padded, 
                         lower=1./maxperiod, 
                         upper=1./minperiod,
                         num= 2000)
            px = np.sqrt(px)
            # look at 20 days or at most at the length of dataset
            pdm_periods, pdm_thetas = pdm(self.t, self.m, minperiod, maxperiod, 0.5/24)
            period = pdm_periods[np.argmin(pdm_thetas)]
            psd_period = 1./f[np.argmax(px)]
            #psd_freq = f[np.argmax(px)]
            
            from scipy import interpolate
            
            i = np.argsort(1./f)
            psd_periods = 1./f[i]
            psd_power = px[i]
            psd_periods = np.insert(psd_periods, 0, 0.0)
            psd_power = np.insert(psd_power,0 ,0.0)
            psd_int = interpolate.interp1d(psd_periods, psd_power)
            
            # use interpolation function returned by `interp1d`
            sum_amp = psd_int(pdm_periods)*(1.-pdm_thetas)   
            i = np.argmax(sum_amp)
            period = pdm_periods[i] 
            theta = pdm_thetas[i]
            
            period = star['period']
            period_err = star['period_err']
            tp, yp = phase(self.t, self.m, period)

                
            s1 = np.sin(2*np.pi*tp/period)
            c1 = np.cos(2*np.pi*tp/period)
            s2 = np.sin(4*np.pi*tp/period)
            c2 = np.cos(4*np.pi*tp/period)
            
            A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
            c, resid,_,_ = np.linalg.lstsq(A,yp)
            amp_err = resid[0]

            amp = max(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)-\
                  min(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)

            tp1 = np.linspace(0.0, period, 100)
            s1 = np.sin(2*np.pi*tp1/period)
            c1 = np.cos(2*np.pi*tp1/period)
            s2 = np.sin(4*np.pi*tp1/period)
            c2 = np.cos(4*np.pi*tp1/period)
                


            plt.subplot(411) ##################################################
            plt.title('%s (%d) B-V=%.2f' % (starid, star['tab'], star['bv']))
            self.plot_lightcurve()
            
            plt.subplot(412) ##################################################
            plt.axvline(x = psd_period, color='green', alpha=0.5)
            plt.axvline(x = period, color='red', alpha=0.5)
            plt.semilogx(1./f,px*1000, 'k')
            plt.xlim(0.1, 30)
            plt.grid()

            plt.subplot(413) ##################################################
            plt.plot(pdm_periods, pdm_thetas, 'k')
            from functions import normalize
            plt.plot(pdm_periods, normalize(sum_amp), 'b')
            
            #plt.ylim(theta, 1.0)
            plt.axvline(x = psd_period, color='green')
            plt.axvline(x = period, color='red')
            plt.grid()
            
            
            plt.subplot(414) ##################################################
            plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
            plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                     linestyle='--', linewidth=2)
            plt.xlim(0.0,period)
            plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
            plt.xlabel('P = %.4f' % period)
            plt.grid()
            #plt.show()
            comment = 'P=%6.3f+-%.3f a=%.3f+-%.4f %.2f' % \
            (period, period_err, amp, amp_err, theta)
            #plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
            pdf.savefig()
            plt.close()
                
            logger.info( comment)
            print comment
            
if __name__ == '__main__':
    m48 =  M48ReAnalysis(config.datapath)
    m48.getstars()
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(config.plotpath+'lightcurves.pdf') as pdf:
        m48.analysis()
        