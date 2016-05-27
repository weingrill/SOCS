#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 27, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
import config

class NGC6633gab(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def loadfromfile(self, filename= '/work2/jwe/Projects/NGC6633/data/tallrefavgcal_fields.dat'):
        self.table = np.genfromtxt(filename, names = True, missing_values='nan')
        i = np.where(self.table['Bmag'] > 25.0)
        self.table['Bmag'][i] = np.nan
        i = np.where(self.table['Vmag'] > 25.0)
        self.table['Vmag'][i] = np.nan
        i = np.where(self.table['Imag'] > 25.0)
        self.table['Imag'][i] = np.nan
        
        maxsigma = 0.05
        i = np.where(self.table['s_Bmag'] < maxsigma)
        self.table = self.table[i]
        i = np.where(self.table['s_Vmag'] < maxsigma)
        self.table = self.table[i]
        i = np.where(self.table['s_Imag'] < maxsigma)
        self.table = self.table[i]
        i = np.where(np.isfinite(self.table['RA_J2000']) & np.isfinite(self.table['DEC_J2000']))
        self.table = self.table[i]
        print self.table.shape
        
    def plot(self):
        import matplotlib.pyplot as plt
        
        ra = self.table['RA_J2000']
        dec = self.table['DEC_J2000']
        mra = np.mean(self.table['RA_J2000'])
        mde = np.mean(self.table['DEC_J2000'])
        
        d = np.sqrt((ra-mra)**2 + (dec-mde)**2)
        print np.nanmin(d), np.nanmax(d)
        bv = self.table['Bmag']-self.table['Vmag']
        vi = self.table['Vmag']-self.table['Imag']
        vmag = self.table['Vmag']
        c = vi
        
        # vary the transparency with distance
        min_dist = 5.0
        max_dist = 20.0
        from_dist = np.linspace(max_dist,min_dist,10)
        delta_dist = from_dist[0] - from_dist[1]
        to_dist = from_dist + delta_dist
        alphas = np.linspace(0.0,1.0, 10)
        
        for i in range(10):
            print from_dist[i],to_dist[i],alphas[i]
            j = np.where((d* 60. > from_dist[i]) & (d* 60. <= to_dist[i]))
            plt.scatter(bv[j], vmag[j],marker='.',edgecolor='none',c=c[j], alpha=alphas[i])
        
        #plot the core of the cluster intransparent
        i = np.where(d* 60. < min_dist)
        print i
        plt.scatter(bv[i], vmag[i],marker='.',edgecolor='none',c=c[i], alpha=1.0)
        
            
        plt.xlim(0.0,2.0)
        plt.ylim(20.0, 12.0)
        plt.savefig(config.resultpath+'ngc6633gab_cmd.pdf')
        plt.close()
    
    def plot_colorcolor(self):
        import matplotlib.pyplot as plt
        from functions import scaleto
        iso = np.genfromtxt('/home/jwe/data/iso_00.500_Gyr.dat', names=True)
        ra = self.table['RA_J2000']
        dec = self.table['DEC_J2000']
        mra = np.mean(self.table['RA_J2000'])
        mde = np.mean(self.table['DEC_J2000'])
        d = np.sqrt((ra-mra)**2 + (dec-mde)**2)*60.
        
        bv = self.table['Bmag']-self.table['Vmag']
        vi = self.table['Vmag']-self.table['Imag']
        vmag = self.table['Vmag']
        s = scaleto(vmag, (40.,1.))
        i = np.where(d < 8.)
        plt.scatter(bv[i]-0.182, vi[i]-0.3,marker='o',edgecolor='none',c=d[i], s=s[i], alpha=0.25)
        #plt.scatter(iso['BV'],iso['VI'],marker='o',edgecolor='none', s=scaleto(iso['M_V'],(40.1,1.)), alpha=0.5)
        plt.plot(iso['BV'],iso['VI'],'k')
        #plt.colorbar()
        plt.xlabel('B - V')
        plt.ylabel('V - I')
        
        plt.savefig(config.resultpath+'ngc6633gab_bv_vi.pdf')
        plt.close()
        
    def savetofile(self, filename= '/work2/jwe/Projects/NGC6633/data/ngc6633gab.dat'):
        np.savetxt(filename, self.table, fmt='%d %.6f %.6f %.3f %.3f %.3f %.3f %.3f %.3f', header='ID RA_J2000 DEC_J2000 Imag s_Imag Vmag s_Vmag Bmag s_Bmag')
        
if __name__ == '__main__':
    n = NGC6633gab()
    n.loadfromfile()
    #n.plot()
    n.plot_colorcolor()
    n.savetofile()
    