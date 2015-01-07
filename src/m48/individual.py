'''
Created on Dec 11, 2014

@author: jwe
'''
from m48 import M48Star
from functions import sigma_clip, phase
import config
import pylab as plt
import numpy as np

class M48Analysis(M48Star):
    '''
    Analysis for individual objects
    '''


    def _load_lightcurve(self):
        t, m, _ = self.lightcurve()
        t -= min(t)
        
        # perform a 3sigma clipping
        self.t, self.m = sigma_clip(t, m)

    
    def _init_plot(self):
            
        from matplotlib import rcParams
        
        fig_width = 18/2.54  # width in inches, was 7.48in
        fig_height = 26/2.54  # height in inches, was 25.5
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

    def plot_lightcurve(self, show=False):
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
#        if show: plt.show()
#        plt.savefig(config.plotpath+'%s_lightcurve.pdf' % self.starid)
#        plt.close()
     
    def plot_psd(self, show=False):
        from psd import ppsd
        from functions import gauss_fit
        # perform a power spectrum analysis
        t, m = self.t, self.m- np.mean(self.m)
        n = len(t)
        t_padded = np.zeros(4*n)
        t_padded[:n] = t
        t_padded[n:] = np.linspace(max(t),4*max(t),3*n)
        m_padded = np.zeros(4*n)
        m_padded[:n] = m
        
        px, f = ppsd(t_padded, m_padded, lower=1./20, upper=1./0.1)
        px = np.sqrt(px)
        period1 = 1./f[np.argmax(px)]
        popt = gauss_fit(f, px, amp=max(px), mean=1.0/period1, sigma=0.1)
        print popt
        t_window = np.linspace(0,4*max(t),4*max(t)*24)
        m_window = np.zeros(len(t_window))
        for tt in t:
            i = np.argmin(abs(t_window-tt))
            m_window[i] = max(abs(m))/2
            
        p_win = abs(np.fft.rfft(m_window))
        p_win /= p_win.size
        f_win = np.fft.fftfreq(t_window.size, 1.0/24)[:p_win.size]
        
        plt.axvline(x = period1, color='green', alpha=0.5)
        #plt.semilogx(1./f,px, 'k')
        plt.plot(1./f,px*1000.0, 'k')
        plt.plot(1./f_win[1:],p_win[1:]*1000.0, 'r')
        #plt.plot(t_window,m_window)
        plt.xlim(0.1, 20)
        plt.xlabel('period [days]')
        plt.ylabel('semi-amplitude [mmag]')
        plt.grid()
#        if show: plt.show()
#        plt.savefig(config.plotpath+'%s_psd.pdf' % self.starid)
#        plt.close()
        

    def plot_pdm(self, show=False):
        from pdm import pdm
        # look at 20 days or at most at the length of dataset
        length = min([max(self.t), 20.0])
        p1, t1 = pdm(self.t, self.m, 0.1, length, 1.0/24)
        period = p1[np.argmin(t1)]
        theta = min(t1)
        plt.plot(p1, t1, 'k')
        plt.ylim(theta, 1.0)
        plt.axvline(x = period, color='red')
        plt.xlabel('period [days]')
        plt.ylabel('theta')
        plt.grid()
#        if show: plt.show()
#        plt.savefig(config.plotpath+'%s_pdm.pdf' % self.starid)
#        plt.close()

    def plot_lsq(self, show=False):
        from analysis import lsqspectrum
        t,l = self.t, self.m-np.mean(self.m)
        lsq_amp,lsq_p, lsq_pha, _, _, res = lsqspectrum(t, l, limit = 30)
        print np.mean(res)
        print lsq_p
        plt.vlines(lsq_p, 0., lsq_amp*1000.0*10)
        #plt.xlim(0.1, 20)
        plt.xlabel('period [days]')
        plt.ylabel('amplitude [mmag]')
        
        x = np.linspace(0.0, 62.0, 200)
        y = np.zeros(200)
        for a,p,pha in zip(lsq_amp,lsq_p, lsq_pha):
            y += a*np.cos(2*np.pi/p*x+pha)
        plt.plot(x,y*1000,'r')
        plt.scatter(t,l*1000,edgecolor='none',facecolor='g')    
#        if show: plt.show()
#        plt.savefig(config.plotpath+'%s_lsq.pdf' % self.starid)
#        plt.close()

    def phase_plot(self, show=False):
        
        period = self.period
        tp, yp = phase(self.t,self.m, period)
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, resid,rank,sigma = np.linalg.lstsq(A,yp)
        print self.starid, c, resid, rank, sigma

        
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
        plt.xlabel('period [days]')
        plt.ylabel('theta')
        plt.grid()
#        if show: plt.show()
#        plt.savefig(config.plotpath+'%s_phase.pdf' % self.starid)
#        plt.close()

if __name__ == '__main__':
    star = M48Analysis('20140303A-0074-0013#1952')
    star._load_lightcurve()
    star._init_plot()
    plt.subplot(411)
    star.plot_lightcurve()
    plt.subplot(412)
    star.plot_psd()
    plt.subplot(413)
    star.plot_pdm()
    plt.subplot(414)
    star.plot_lsq()
    #star.phase_plot()        
    plt.savefig(config.plotpath+'%s.pdf' % star.starid)
    plt.close()