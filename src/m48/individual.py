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
        fig_height = 28/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 7,
          'axes.titlesize': 8,
          'font.size': 7,
          'xtick.labelsize': 7,
          'ytick.labelsize': 7,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)

    def plot_lightcurve(self, show=False, save = False):
        """
        plot the lightcurve for a given star
        """

        mean = np.mean(self.m)
        std = np.std(self.m)
        
        plt.hlines(mean,min(self.t),max(self.t),linestyle='--')
        plt.ylim(mean+std*3,mean-std*3)
        plt.xlim(min(self.t),max(self.t))
        plt.grid()
        plt.title(self.starid+' = star #'+str(self['tab']))
        plt.xlabel('days')
        plt.ylabel('mag')
        plt.scatter(self.t, self.m, edgecolor='none')
        plt.plot(self.t, self.m, 'gray')
        plt.minorticks_on()
     
    def plot_psd(self, show=False, save=False):
        from psd import ppsd
        
        # perform a power spectrum analysis
        t, m = self.t, self.m- np.mean(self.m)
        n = len(t)
        t_padded = np.zeros(8*n)
        t_padded[:n] = t
        t_padded[n:] = np.linspace(max(t),8*max(t),7*n)
        m_padded = np.zeros(8*n)
        m_padded[:n] = m
        
        px, f = ppsd(t_padded, m_padded, lower=1./20, upper=1./0.1)
        px = np.sqrt(px)
        i = np.argmax(px)
        period1 = 1./f[i]
        
        t_window = np.linspace(0,8*max(t),8*max(t)*24)
        m_window = np.zeros(len(t_window))
        for tt in t:
            i = np.argmin(abs(t_window-tt))
            m_window[i] = max(abs(m))/2
            
        p_win = abs(np.fft.rfft(m_window))
        p_win /= p_win.size
        f_win = np.fft.fftfreq(t_window.size, 1.0/24)[:p_win.size]
        
        plt.axvline(x = period1, color='green', alpha=0.5)
        plt.plot(1./f[1:],px[1:]*1000.0, 'r')
        plt.plot(1./f_win[1:],p_win[1:]*1000.0, 'k')
        
        plt.axhline(np.std(px)*5000.0, linestyle='--', color='b')
        
        plt.xticks(np.arange(20))
        plt.xlim(0.1, 20)
        plt.xlabel('period [days]')
        plt.ylabel('semi-amplitude [mmag]')
        plt.minorticks_on()
        plt.grid()
        if show: plt.show()
        elif save: 
            plt.savefig(config.plotpath+'%s_psd.pdf' % self.starid)
        
    def plot_pdm(self, show=False, save = False):
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
        plt.minorticks_on()
        plt.grid()

    def plot_lombscargle(self):
        from scipy.signal import lombscargle
        
        n = 2000
        f = np.linspace(2.*np.pi/20, 2.*np.pi/0.1, n)
        pgram = lombscargle(self.t, self.m, f)
        period = 1./(f/(2.*np.pi))
        plt.plot(period, np.sqrt(4*(pgram/n)), 'r')
        plt.xlabel('period [days]')
        plt.ylabel('FAP')
        plt.xticks(np.arange(20))
        plt.xlim(0.1, 20)
        plt.minorticks_on()
        plt.grid()
        

    def plot_lsq(self, show=False, save = False):
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

    def phase_plot(self, show=False, save = False):
        
        period = self['p_fin']
        tp, yp = phase(self.t,self.m, period)
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, _,_,_ = np.linalg.lstsq(A,yp)
                
        plt.scatter(tp, yp-np.mean(yp), edgecolor='none', alpha=0.75)
        plt.scatter(tp+period, yp-np.mean(yp), edgecolor='none', alpha=0.75)
        tp1 = np.linspace(0.0, period*2.0, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)
        
        plt.axvline(period, linestyle='-.')
        plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                 linestyle='--', linewidth=2)
        plt.plot(tp1,c[1]*s1, 'r', linewidth=0.5)
        plt.plot(tp1,c[2]*c1, 'g', linewidth=0.5)
        plt.plot(tp1,c[3]*s2, 'y', linewidth=0.5)
        plt.plot(tp1,c[4]*c2, 'm', linewidth=0.5)
        plt.xlim(0.0,period*2)
        plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
        plt.xlabel('period [days]')
        plt.ylabel('mag')
        plt.minorticks_on()
        plt.grid()

    def lightcurve_overplot(self, show=False, save = False):
        
        period = self['period']
        tp, yp = phase(self.t,self.m, period)
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, resid,rank,sigma = np.linalg.lstsq(A,yp)
        print self.starid, c, resid, rank, sigma

        
        plt.scatter(self.t,self.m-np.mean(self.m), edgecolor='none', alpha=0.75, c='red')
        tp1 = np.linspace(self.t[0], self.t[-1], 200)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)
        
        
        plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', linewidth=1)
        
        plt.xlim(0.0,tp1[-1])
        plt.ylim(0.03,-0.03)
        plt.title('star #%d' % self['tab'])
        plt.xlabel('time [days]')
        plt.ylabel('V [mag] + %.2f' % np.mean(self.m))
        plt.grid()
        if show: plt.show()
        elif save: 
            plt.savefig(config.plotpath+'%s_lcsine.pdf' % self.starid)
            plt.savefig(config.plotpath+'%s_lcsine.eps' % self.starid)

    def plot_clean(self, show=False, save = False):
        a = np.loadtxt(config.datapath+'/clean/%s.ncfile' % self.starid[7:])
        f = a[:,0]
        px = a[:,1]*1000.0
        
        #px = np.sqrt(px*1000.0)/2.0
        
        plt.plot(1./f[1:],px[1:], 'r')
        plt.xticks(np.arange(20))
        plt.xlim(0.1, 20)
        plt.minorticks_on()
        plt.xlabel('period [days]')
        plt.ylabel('power')
        plt.axvline(1./f[np.argmax(px)], color='green', alpha=0.5)
        plt.axhline(np.std(px)*5.0, color='blue', linestyle='--')
        plt.grid()
        if show: plt.show()
        elif save: 
            plt.savefig(config.plotpath+'%s_clean.pdf' % self.starid)
            plt.savefig(config.plotpath+'%s_clean.eps' % self.starid)
        

if __name__ == '__main__':
    from matplotlib.backends.backend_pdf import PdfPages
    starids = """20140303A-0074-0013#1952
20140306A-0097-0013#199
20140307A-0000-0013#787
20140301A-0071-0013#1178
20140302A-0000-0015#597
20140303A-0074-0013#1855
20140301A-0071-0013#1316
20140301A-0071-0013#1218
20140303A-0074-0013#1062
20140302A-0000-0015#1449
20140307A-0000-0013#342
20140306A-0097-0013#712
20140302A-0000-0015#1091
20140306A-0097-0013#239
20140303A-0074-0013#1041
20140306A-0097-0013#988
20140307A-0000-0013#39
20140306A-0097-0013#880
20140301A-0071-0013#125
20140307A-0000-0013#1403
20140303A-0074-0013#1063
20140306A-0097-0013#802
20140302A-0000-0015#663
20140307A-0000-0013#284
20140307A-0000-0013#486
20140306A-0097-0013#329
20140306A-0097-0013#115
20140306A-0097-0013#703
20140302A-0000-0015#928
20140306A-0097-0013#256
20140307A-0000-0013#436
20140307A-0000-0013#202
20140301A-0071-0013#1111
20140303A-0074-0013#1745
20140306A-0097-0013#546
20140307A-0000-0013#1466
20140306A-0097-0013#716
20140306A-0097-0013#1796
20140306A-0097-0013#1022
20140306A-0097-0013#1807
20140302A-0000-0015#524
20140307A-0000-0013#1529
20140307A-0000-0013#334
20140306A-0097-0013#732
20140303A-0074-0013#1850
20140307A-0000-0013#735
20140303A-0074-0013#1414
20140301A-0071-0013#1349
20140303A-0074-0013#1942
20140301A-0071-0013#329
20140302A-0000-0015#522
20140306A-0097-0013#525
20140306A-0097-0013#1167
20140307A-0000-0013#684
20140306A-0097-0013#376
20140306A-0097-0013#660
20140307A-0000-0013#220
20140306A-0097-0013#357
20140302A-0000-0015#6
20140306A-0097-0013#1359
20140306A-0097-0013#406
20140306A-0097-0013#1336
20140302A-0000-0015#176
20140306A-0097-0013#745
20140301A-0071-0013#357
20140301A-0071-0013#1614
20140306A-0097-0013#1823
20140306A-0097-0013#270
20140307A-0000-0013#1514
20140302A-0000-0015#488"""

    pdf = PdfPages(config.plotpath+'M48analysis.pdf')
    for starid in starids.split('\n')[1:5]:
        print starid
        star = M48Analysis(starid)
        star._load_lightcurve()
        star._init_plot()
        plt.subplot(511)
        star.plot_lightcurve()
        plt.subplot(512)
        star.plot_lombscargle()
        plt.subplot(513)
        star.plot_clean()
        plt.subplot(514)
        star.plot_pdm()
        plt.subplot(515)
        
        star.phase_plot()        
        pdf.savefig()
        plt.close()
        #plt.savefig(config.plotpath+'%s.pdf' % star.starid)
    pdf.close()