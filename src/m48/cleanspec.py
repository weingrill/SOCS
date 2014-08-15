'''
Created on Aug 8, 2014

@author: jwe
'''
from m48star import M48Star


def loadspectra():
    path = '/work1/jwe/Dropbox/M48/data/spectra_syd'
    
    from glob import glob
    import numpy as np
    import os.path
    #from datasource import DataSource
    import matplotlib.pyplot as plt
    from functions import gauss_fit, gauss
    #wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
    
    files = glob(path+'/*.ncfile')
    flen = len(files)
    bvarr = np.zeros(flen)
    freqarr = np.zeros(flen)
    
    lf = open('/work1/jwe/Dropbox/M48/data/cleanspec.dat','wt')
    
    for f in files[:flen]:
        shortid = '%' + os.path.basename(f).rstrip('.ncfile')
        star = M48Star(shortid)
        i = files.index(f)
        print star.starid,
        a = np.loadtxt(f)
        freq = a[:1600,0]
        amp = a[:1600,1]
        
        m = np.argmax(amp)
        n = len(freq)
        assert(n<=1600)
        kleft=m
        while amp[kleft-1]<amp[kleft] and kleft>0: kleft -= 1
        kright=m
        while kright+1<n and amp[kright+1]<amp[kright]: kright += 1
        
        sfreq = freq[kleft:kright]
        samp = amp[kleft:kright]
        print '%.3f' % star['bv'],
        
        ma = np.mean(amp)
        speriod = 1./sfreq
        j = np.argsort(speriod)
        speriod = speriod[j]
        samp = samp[j]
    
        mean_period = 1./freq[m]
        sig = (max(speriod)-min(speriod))/2.0 #np.sqrt(mean_period)
        
        g = gauss_fit(speriod, samp-ma, amp=amp[m]-ma, mean=mean_period, sigma = sig)
        freqarr[i] = g[1]
        print '%.5f %.5f %.5f' % (g[1], g[0], g[2])
        bvarr[i] = star['bv']
        lf.write('%25s %.3f %.5f %.5f %.5f\n' % (star.starid, star['bv'], g[1], g[0], g[2]))
        
        #plt.plot(speriod,samp,'b')
        #plt.plot(speriod,gauss(speriod,*g)+ma,'r',label='fit')
        #plt.show()
        
    lf.close()
    plt.scatter(bvarr, 1./freqarr)
    plt.grid()
    plt.xlim(0.4, 2.0)
    plt.ylim(0.5, 30)
    plt.xlabel('B-V')
    plt.ylabel('CLEAN period (days)')
    plt.show()

def inject_periods():
    
    lf = open('/work1/jwe/Dropbox/M48/data/cleanspec.dat','rt')
    lines = lf.readlines()
    lf.close()
    for l in lines:
        ls = l.split()
        starid = ls[0].strip()
        star = M48Star(starid)
        print starid,ls[2:5]
        star['clean_period'] = ls[2]
        star['clean_amp'] = ls[3]
        if not ls[4]=='inf':
            star['clean_sigma'] = ls[4]
        
if __name__ == '__main__':
    inject_periods()