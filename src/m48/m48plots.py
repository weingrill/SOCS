'''
Created on Sep 18, 2014

@author: jwe
'''
from m48star import M48Star            
import numpy as np
import pylab as plt
from pdm import pdm
from psd import ppsd
from matplotlib import rcParams
from functions import sigma_clip, phase

class M48Plots(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.stars = []
        self.getstars()

    def getstars(self):
        """
        build up a list of stars, where we do not have periods yet
        """
        
        query = """SELECT starid 
        FROM m48stars 
        WHERE pman>0;"""
        
        result = self.wifsip.query(query)
        print '... %d stars found' % len(result)
        self.stars = [s[0] for s in result]

    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        
        mean = np.mean(self.m)
        #std = np.std(self.m)
        plt.hlines(mean,min(self.t),max(self.t),linestyle='--')
        #plt.ylim(mean+std*3,mean-std*3)
        plt.xlim(min(self.t),max(self.t))
        plt.grid()
        #plt.scatter(self.t, self.m, edgecolor='none')
        plt.errorbar(self.t, self.m, yerr=self.e*0.5, fmt='o')
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])

    def phase_plot(self, starid):
        star = M48Star(starid)
        try:
            t, m, e = star.lightcurve()
        except TypeError:
            print 'no data'
            return
        t, m, e = sigma_clip(t, m, e)
        fperiod = star['period']
        clean_period = star['clean_period']
        pman = star['pman']
        if abs(fperiod-pman)/pman < abs(clean_period-pman)/pman:
            period = fperiod
        else:
            period = clean_period
        if abs(period-pman)>1:
            period = pman
        print '%.3f (%.2f)' % (period, pman)
        tp, yp = phase(t, m, period)
        plt.scatter(tp, yp-np.mean(yp), edgecolor='none', facecolor='k', size=5)
        plt.scatter(tp+period, yp-np.mean(yp), edgecolor='none', facecolor='k', size=5)
        plt.axvline(period, linestyle='--', color='k')
        plt.title(starid)
        plt.xlabel('period %.3f (%.2f) days' % (period, pman))
        plt.ylabel(star['vmag'])
        plt.xlim(0,period*2)
        plt.ylim(plt.ylim()[1],plt.ylim()[0])
        plt.grid()
        
        #plt.close()

    def plot(self, show=False):
        """plot lightcurve"""
        
        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 10,
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

        
        sp = 1
        phase=1
        for starid in self.stars:
            print starid,
            plt.subplot(4,2,sp)
            sp += 1
            self.phase_plot(starid)
            if sp==9 or starid==self.stars[-1]:
                plt.tight_layout()
                #plt.show()
                plt.savefig('/work2/jwe/m48/plots/phase%d.pdf' % phase)
                phase += 1
                sp = 1 
                plt.close()


m48plots = M48Plots()
m48plots.plot(show=True)