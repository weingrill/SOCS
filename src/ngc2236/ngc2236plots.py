#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
import numpy as np
from matplotlib import pyplot as plt

class Ngc2236Plots(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        wifsip = DataSource(database=config.dbname, 
                                 user=config.dbuser, 
                                 host=config.dbhost,
                                 dictcursor=True)
        query = """SELECT bmag, vmag, period, good, member 
        FROM ngc2236
        WHERE NOT bv is NULL AND NOT vmag is NULL"""
        self.data = wifsip.query(query) 
        self.dm = 14.0
         

    def plot_cpd(self):
        vmag = np.array([d['vmag'] for d in self.data])
        bmag = np.array([d['bmag'] for d in self.data])
        period = np.array([d['period'] for d in self.data])
        good = np.array([d['good'] for d in self.data])
        member = np.array([d['member'] for d in self.data])
        vmag = 1.0001722499859311*vmag + 0.5402413960734754
        bmag = 1.0108861546902421*bmag + 0.9271685986891227
        bv = bmag-vmag
        
        i = np.where(period>0.0)
        plt.scatter(bv[i]-0.56, period[i], facecolor='none', edgecolor='k',label='WiFSIP')
        j = np.where((period>0.0) & (good==True))
        plt.scatter(bv[j]-0.56, period[j], facecolor='r', edgecolor='none',label='rotators')
        k = np.where(member==True)
        plt.scatter(bv[k]-0.56, period[k], marker='s', s=50, facecolor='none', edgecolor='g',label='members')
        plt.ylim(0.0,15.0)
        plt.xlim(-0.5, 2.0)
        plt.xlabel('(B - V)$_0$')
        plt.ylabel('P [days]')

        import gyroage
        
        bvgyro = np.linspace(0.4, 1.5, num=50)
        #bv = np.linspace(0.473, 1.631, num=100)
        
        P500 = gyroage.gyroperiod(bvgyro, 500)
        plt.plot(bvgyro, P500, color='r')
        
        P800 = gyroage.gyroperiod(bvgyro, 800)
        plt.plot(bvgyro, P800, 'b--')
        
        plt.savefig(config.plotpath+'NGC2236cpd.eps')
        plt.show()
        
    
    def plot_cmd(self):
        vmag = np.array([d['vmag'] for d in self.data])
        bmag = np.array([d['bmag'] for d in self.data])
        period = np.array([d['period'] for d in self.data])
        good = np.array([d['good'] for d in self.data])
        member = np.array([d['member'] for d in self.data])
        # calibration using data from Babu
        vmag = 1.0001722499859311*vmag + 0.5402413960734754
        bmag = 1.0108861546902421*bmag + 0.9271685986891227
        bv = bmag-vmag
        
        plt.scatter(bv, vmag, facecolor='none', edgecolor='k',label='WiFSIP')
        i = np.where((period>0.0) & (good==True))
        plt.scatter(bv[i], vmag[i], facecolor='r', edgecolor='none',label='rotators')
        j = np.where(member==True)
        plt.scatter(bv[j], vmag[j], marker='s', s=50, facecolor='none', edgecolor='g',label='members')

        iso_v, iso_bv = self.load_isochrone() 
        plt.plot(iso_bv+0.56, iso_v, 'b', alpha=0.3, lw=5.0,label='500 Myr iso')
        
        plt.ylim(17.0, 12.0)
        plt.xlim(0.0, 2.0)
        plt.xlabel('B - V')
        plt.ylabel('V mag')
        plt.legend()
        plt.savefig(config.plotpath+'NGC2236cmd.eps')
        plt.show()

    def load_isochrone(self):
        from numpy import loadtxt
        isofile = config.datapath+'0p500Gyr_FeH0p0_Y0p277_AMLTsol.iso'
        a = loadtxt(isofile)
        iso_mv = a[:,5]
        iso_bv = a[:,6]
        return iso_mv+self.dm, iso_bv
     
if __name__ == '__main__':
    nplt = Ngc2236Plots()
    nplt.plot_cmd()
    nplt.plot_cpd()