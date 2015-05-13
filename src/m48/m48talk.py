'''
Created on Jul 14, 2014

@author: jwe
'''

from m48 import M48Star
from functions import sigma_clip, phase
import numpy as np
import pylab as plt


class M48TalkStar(M48Star):
    def _load_lightcurve(self):
        t, m, _ = self.lightcurve()
        t -= min(t)
        
        # perform a 3sigma clipping
        self.t, self.m = sigma_clip(t, m)

    
    def _init_plot(self):
            
        from matplotlib import rcParams
        
        fig_width = 12/2.54  # width in inches, was 7.48in
        fig_height = 8/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)

    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        mean = np.mean(self.m)
        std = np.std(self.m)
        
        plt.hlines(mean,min(self.t),max(self.t),linestyle='--')
        plt.ylim(mean+std*3,mean-std*3)
        plt.xlim(min(self.t),max(self.t))
        plt.grid()
        plt.title(self.starid)
        plt.xlabel('days')
        plt.ylabel('mag')
        plt.scatter(self.t, self.m, edgecolor='none')
        #plt.show()
        plt.savefig('/work1/jwe/Dropbox/Documents/Talks/AIP2014a/lightcurve.pdf')
        plt.close()
     
    def plot_psd(self):
        from psd import ppsd
        # perform a power spectrum analysis

        px, f = ppsd(self.t, self.m- np.mean(self.m), lower=1./30, upper=1./0.1)
        px = np.sqrt(px)
        period1 = 1./f[np.argmax(px)]
        
        
        plt.axvline(x = period1, color='green', alpha=0.5)
        plt.semilogx(1./f,px, 'k')
        plt.xlim(0.1, 30)
        plt.title(self.starid)
        plt.xlabel('period [days]')
        plt.ylabel('semi-amplitude mag')
        plt.grid()
        plt.show()
        plt.savefig('/work1/jwe/Dropbox/Documents/Talks/AIP2014a/psd.pdf')
        plt.close()
        

    def plot_pdm(self):
        from pdm import pdm

        # look at 20 days or at most at the length of dataset
        length = min([max(self.t), 20.0])
        p1, t1 = pdm(self.t, self.m, 0.1, length, 60./86400.)
        period = p1[np.argmin(t1)]
        theta = min(t1)
        plt.plot(p1, t1, 'k')
        plt.ylim(theta, 1.0)
        plt.axvline(x = period, color='red')
        plt.title(self.starid)
        plt.xlabel('period [days]')
        plt.ylabel('theta')
        plt.grid()
        plt.show()
        plt.savefig('/work1/jwe/Dropbox/Documents/Talks/AIP2014a/pdm.pdf')
        plt.close()

    def phase_plot(self):
        period = self.period
        tp, yp = phase(self.t,self.m, period)
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, _,_,_ = np.linalg.lstsq(A,yp)
        
        plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
        tp1 = np.linspace(0.0, period, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)
        
        
        plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                 linestyle='--', linewidth=2)
        plt.plot(tp1,c[1]*s1, 'r', linewidth=0.5)
        plt.plot(tp1,c[2]*c1, 'g', linewidth=0.5)
        plt.plot(tp1,c[3]*s2, 'y', linewidth=0.5)
        plt.plot(tp1,c[4]*c2, 'm', linewidth=0.5)
        plt.xlim(0.0,period)
        plt.ylim(0.05,-0.05)
        plt.title(self.starid)
        plt.xlabel('period [days]')
        plt.ylabel('theta')
        plt.grid()
        #plt.show()
        plt.savefig('/work1/jwe/Dropbox/Documents/Talks/AIP2014a/phase.pdf')
        plt.close()
        
    def plot_cpd(self):
        '''
        create plot for Heraeus Talk
        '''
        from datasource import DataSource
        wifsip = DataSource(database='wifsip', user='sro', host='oldpina.aip.de')
        query = """SELECT bv, period 
                    FROM m48stars 
                    WHERE good;"""

        data = wifsip.query(query)
        
        bv = np.array([d[0] for d in data])
        period = np.array([d[1] for d in data])
        
        import gyroage
        from functions import logspace
        
        bv360 = logspace(0.5, 2.0, num=100)
        #P = gyroage.gyroperiod(bv360, 360.0, version=2007)
        P, pc = gyroage.gyroperiod(bv360, 360.0, version=2003)
        plt.plot(bv360, pc, color='b', linestyle='--')
        plt.plot(bv360, P, color='r')
        
        plt.scatter(bv-self.ebv, period, s=1, 
                    edgecolor='none', c='k')
         
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('period [days]')
        plt.ylim(0.0, 20.0)
        plt.xlim(0.0, 2.0)
        plt.grid()
        plt.savefig('/home/jwe/Documents/Talks/Heraeus2015/m48cpd.eps')
        plt.savefig('/home/jwe/Documents/Talks/Heraeus2015/m48cpd.pdf')
        plt.close()

        


if __name__ == '__main__':
    star = M48TalkStar('20140302A-0000-0015#278')
    star.plot_cpd()
    #star._load_lightcurve()
    #star._init_plot()
    
    #star.plot_lightcurve()
    #star.plot_psd()
    star.phase_plot()