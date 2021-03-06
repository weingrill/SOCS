'''
Created on Aug 8, 2014

@author: jwe
'''
from m48star import M48Star
import numpy as np

def loadspectra():
    path = '/work1/jwe/Dropbox/M48/data/spectra_syd'
    
    from glob import glob
    import numpy as np
    import os.path
    #from datasource import DataSource
    import matplotlib.pyplot as plt
    from functions import gauss_fit
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

def update_periodic():
    lf = open('/work2/jwe/m48/data/list_periodic_ms.txt','rt')
    lines = lf.readlines()
    lf.close()
    for l in lines:
        if not l[0]=='#':
            ls = l.split('\t')
            star = M48Star('%'+ls[0])
            manualP = float(ls[1])
            quality = ls[2].rstrip(' \n')
            if len(ls)>=4:
                notes = ls[3].rstrip(' \n')
            else:
                notes = None 
            print '%-25s\t%.3f\t%.2f\t%s\t%s' % (star['starid'], 
                        star['clean_period'], 
                        manualP, 
                        quality, 
                        notes)
            star['pman'] = manualP
            star['quality'] = quality.rstrip()
            star['notes'] = notes
   
def plot_clean():
    import pylab as plt
    import numpy as np
    from datasource import DataSource
    wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
    query = "SELECT vmag+0.06, bv FROM m48stars WHERE NOT bv is NULL;"
    data = wifsip.query(query)
    vmag = np.array([d[0] for d in data])
    bv = np.array([d[1] for d in data])
    
    query = """SELECT vmag+0.06, bv 
                FROM m48stars 
                WHERE pman>0;"""
    data = wifsip.query(query)
    vmag_pman = np.array([d[0] for d in data])
    bv_pman = np.array([d[1] for d in data])
    
    plt.subplot(2,1,1)
    
    plt.scatter(bv, vmag, marker='.',edgecolor='none', facecolor='gray')
    plt.scatter(bv_pman,vmag_pman,marker='o',edgecolor='none', facecolor='g')

    query = """SELECT vmag+0.06, bv, starid 
                FROM m48stars 
                WHERE good;"""
    data = wifsip.query(query)
    vmag_good = np.array([d[0] for d in data])
    bv_good = np.array([d[1] for d in data])
    text = [d[2].split('#')[1] for d in data]
    plt.scatter(bv_good,vmag_good,marker='o',edgecolor='r', facecolor='none')
    for x,y,t in zip(bv_good,vmag_good,text):
        plt.text(x,y,' '+t)
    plt.ylim(16,11.5)
    plt.xlim(0.4,1.15)
    plt.grid()
    
    plt.subplot(2,1,2)
    

    query = """SELECT bv, clean_period, pman 
                FROM m48stars 
                WHERE pman>0;"""
    data = wifsip.query(query)
    bv_syd = np.array([d[0] for d in data])
    clean_syd = np.array([d[1] for d in data])
    pman_syd = np.array([d[2] for d in data])
    plt.scatter(bv_syd,pman_syd,marker='o',edgecolor='none', facecolor='g')
    plt.vlines(bv_syd, clean_syd, pman_syd)
    plt.axhline(1.0, linestyle='-.')

    query = """SELECT bv, clean_period, period, starid 
                FROM m48stars 
                WHERE good;"""
    data = wifsip.query(query)
    bv_good = np.array([d[0] for d in data])
    clean_good = np.array([d[1] for d in data])
    period_good = np.array([d[2] for d in data])
    text = [d[3].split('#')[1] for d in data]
    plt.scatter(bv_good,clean_good,marker='o',edgecolor='r', facecolor='none')
    plt.vlines(bv_good, clean_good, period_good, color='r')
    for x,y,t in zip(bv_good,clean_good,text):
        plt.text(x,y,' '+t)
    plt.ylim(0,10)
    plt.xlim(0.4,1.15)
    plt.grid()
    plt.show()

def redoclean():
    #from datasource import DataSource
    #wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
    #query = """SELECT tab FROM m48stars WHERE good;"""
    #result = wifsip.query(query)
    from starids import starids
    from matplotlib import pyplot as plt
    from functions import gauss, gauss_fit

    for starid in starids.split('\n'):
        doplot = False
        star = M48Star(starid)
        f, a = star.cleanspectrum()
        p = 1./f
        i = np.where((p>1.3) & (p<20.0))
        pi = p[i]
        fi = f[i]
        ai = a[i]
        j = np.argmax(ai)
        pnew = pi[j]
        anew = ai[j]
        pman = star['pman']
        diff = abs(pnew-pman)
        if diff > 0.5:
            k = np.where((p>pman-1.0) & (p<pman+1.0))
            pk = p[k]
            fk = f[k]
            ak = a[k]
            j = np.argmax(ak)
            pnew = pk[j]
            anew = ak[j]
            diff = abs(pnew-pman)
            doplot = True
        
        # find local minima
        h = np.argwhere(pi == pnew)[0][0] 
        h0 = h
        while ai[h0 - 1] <= ai[h0] and h0>0:
            h0 -= 1
        h1 = h
        while ai[h1 + 1] <= ai[h1] and h1<len(ai):
            h1 += 1
        # perform gauss fit
        par = gauss_fit(pi[h0:h1], ai[h0:h1]-min(ai[h0:h1]), amp=anew, mean=pnew, sigma=0.5 )
        #print p
        pgauss = par[1]
        perr = par[2]
        if perr>1.0:
            doplot = True
        # update database
        x = np.linspace(pi[0], pi[-1], 500)
        #gauss(x, a, x0, sigma):
        diff = abs(pnew-pman)
            
        print '%25s %4d %6.3f %5.2f %6.3f %6.3f' % (starid, star['tab'], pman, diff, pgauss, perr)
        star['clean_period'] = pgauss
        star['clean_sigma'] = perr
        if doplot:
            plt.plot(pi, ai)
            plt.axvline(pnew, linestyle='--')
            plt.plot(x, gauss(x, par[0], par[1], par[2]), 'r')
            plt.title('#%d' % star['tab'])
            plt.xlabel('period = %.3f days' % pnew)
            plt.show()
            

    
if __name__ == '__main__':
    #inject_periods()
    #update_periodic()
    #plot_clean()
    redoclean()