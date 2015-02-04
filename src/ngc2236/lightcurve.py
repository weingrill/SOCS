'''
Created on Jan 14, 2013

@author: jweingrill@aip.de
'''
class LightCurve(object):
    """lightcurve class to read CoRoT lightcurves"""

    def __init__(self, filename, channel=None):
        """
        Constructor
        """
        from numpy import diff, where, nan
        import pyfits

        hdulist = pyfits.open(filename)
        hdulist.verify('silentfix')
        self.hdr = hdulist[0].header
        tbdata = hdulist[1].data
        self.time = tbdata.field('DATEHEL')
        status = tbdata.field('STATUS')
        if channel is None:
            if 'WHITEFLUX_IMAG' in tbdata.names:
                self.flux = tbdata.field('WHITEFLUX_IMAG')
            elif 'WHITEFLUX' in tbdata.names:
                self.flux = tbdata.field('WHITEFLUX')
        if channel in ['r','red']:
            if 'REDFLUX' in tbdata.names:
                self.flux = tbdata.field('REDFLUX')
        if channel in ['g','green']:
            if 'GREENFLUX' in tbdata.names:
                self.flux = tbdata.field('GREENFLUX')
        if channel in ['b','blue']:
            if 'BLUEFLUX' in tbdata.names:
                self.flux = tbdata.field('BLUEFLUX')
        
        hdulist.close()
        self.corotid = self.hdr['COROTID']

        i = where(status == 0)
        self.time = self.time[i] + 2451545.0 # 1.1.2000 12:00
        self.flux = self.flux[i]

#        self.err = zeros(len(self.time))
#        for i = 0 to len(self.time):
#            self.err[i] = 
        self.err  = abs(diff(self.flux)/diff(self.time))
        self.winid = self.hdr['RUN_CODE']+'_'+\
                     self.hdr['HLFCCDID'][0:2]+\
                     '_%04d' % self.hdr['WIN_ID']
        self.magr = 0.0
        
        if 'MAGNIT_R' in self.hdr:
            self.magr = self.hdr['MAGNIT_R']
        if 'MAGNIT_V' in self.hdr:
            self.magv = self.hdr['MAGNIT_V']
        else: self.magv = nan    
        
        if 'MAGNIT_B' in self.hdr:
            self.magb = self.hdr['MAGNIT_B']
        else: self.magb = nan
        if self.magb == 'NAN': self.magb = nan
        self.bv = float(self.magb)-float(self.magv)
        #print '%s: %.2f' % (filename, self.bv)
    
    def __len__(self):
        return len(self.time)
        
    def mag(self, zero = 0.0):
        """convert flux to magnitutdes with an optional zero magnitude"""
        from numpy import log10, mean
        
        return -2.512*log10(self.flux/mean(self.flux)) + zero
        
    def normalize(self):
        """normalizes the lightcurve to 1.0; returns the mean"""
        from numpy import mean
        m = mean(self.flux)
        self.flux /= m
        self.err /= m
        return m
    
    def detrend(self, times=None, order=2):
        """
        removes a linear or polynomial trend with the order <order>
        returns the fit
        """
        from numpy import polyfit, poly1d
        #numpy.polynomial.polynomial.polyfit(x, y, deg, rcond=None, full=False, w=weights)
        if times is None:
            times = [self.time[0],self.time[-1]]
        else:
            times.insert(0,0.)
            times.append(self.time[-1])
        t0 = self.time[0]
        for i in range(len(times)-1):
            i0 = self.time.searchsorted(times[i])
            i1 = self.time.searchsorted(times[i+1],'right')
            try:
                pf = polyfit(self.time[i0:i1]-t0, self.flux[i0:i1], order)
            except TypeError:
                print i0,i1
            else:
                p = poly1d(pf)
                self.flux[i0:i1] = self.flux[i0:i1] - p(self.time[i0:i1]-t0) + 1.0
        
    def rebin(self, interval = 512./86400., medianbins=False):
        """
        rebin to new interval using the mean of each bin
        interval determines the length of each bin
        medianbins calculates the median in each bin otherwise the mean is taken
        """
        from numpy import seterr, zeros, isnan, compress, arange, median, mean, std
        data = self.time
        # ...+interval so that the last bin includes the last epoch
        bins = arange(self.time[0], self.time[-1]+interval, interval)
        nbins = len(bins)-1
        t = zeros(nbins)
        f = zeros(nbins)
        e = zeros(nbins)
        # adopted from Ian's Astro-Python Code v0.3
        # http://www.mpia-hd.mpg.de/homes/ianc/python/_modules/tools.html
        # def errxy()
        idx = [[data.searchsorted(bins[i]), \
                data.searchsorted(bins[i+1])] for i in range(nbins)]
        seterr(invalid='ignore')
        if medianbins:
            for i in range(nbins):
                f[i] = median(self.flux[idx[i][0]:idx[i][1]])
                t[i] = median(self.time[idx[i][0]:idx[i][1]])
                e[i] = std(self.flux[idx[i][0]:idx[i][1]])
        else:    
            for i in range(nbins):
                f[i] = mean(self.flux[idx[i][0]:idx[i][1]])
                t[i] = mean(self.time[idx[i][0]:idx[i][1]])
                e[i] = std(self.flux[idx[i][0]:idx[i][1]])
        seterr(invalid='warn')
        valid = ~isnan(t)
        self.flux = compress(valid,f)
        self.time = compress(valid,t)
        self.err = compress(valid,e)
        
    def remedian(self, interval = 512./86400.):
        """
        rebin to new interval using the median of each bin
        """
        self.rebin(interval, medianbins=True)
        
    def interpolate(self, interval = 512./86400.):
        """
        interpolate missing datapoints
        """
        from numpy import arange
        from scipy.interpolate import interp1d
        
        f = interp1d(self.time, self.flux)
        self.time = arange(self.time[0], self.time[-1], interval)
        self.flux = f(self.time)

    def slidingaverage(self, width=11):
        from numpy import ones, mean, empty,hstack
        leftpad = ones(width/2)*self.flux[0]
        rightpad = ones(width-width/2)*self.flux[-1]
        padded = hstack((leftpad, self.flux, rightpad))
        n = len(self.flux)
        result = empty(n)
        for i in range(n):
            result[i] = mean(padded[i:i+width])
        self.flux -= result-1.0
        return result

    def slidingmedian(self, width=11):
        from numpy import ones, median, empty,hstack
        leftpad = ones(width/2)*self.flux[0]
        rightpad = ones(width-width/2)*self.flux[-1]
        padded = hstack((leftpad, self.flux, rightpad))
        n = len(self.flux)
        result = empty(n)
        for i in range(n):
            result[i] = median(padded[i:i+width])
        self.flux -= result-1.0
        return result
    
    def phasefold(self, period):
        t = (self.time-self.time[0]) % period
        i = t.argsort()
        self.time = t[i]
        self.flux = self.flux[i]
        
    def fourier(self, sample= 16):
        """
        calculate Fourier Transform and return frequencies, amplitudes and power
        """
        from math import atan2
        from numpy import mean, fft, hamming
        
        m = mean(self.flux)
        y = self.flux - m
        n = len(y)
        h = hamming(n)
        n2 = sample*n
        ft = fft.rfft(y*h, n2)
        # magic number 3.70, maybe coming from hamming window
        amp = 3.70*abs(ft)/n
        pwr = amp**2.0 # response = 20 * np.log10(mag) # in [dB]
        freq = fft.fftfreq(n2, d=self.dt())
        phase = [atan2(c.real,c.imag) for c in ft]  
        #i = where(freq >= 0.0)
        return freq[:n], amp[:n], pwr[:n], phase[:n]

    def autocorr(self, padding=True):
        """
        calculate the Fourier autocorrelation
        """
        from numpy import mean, fft, real
        m = mean(self.flux)
        n = len(self.flux)
        k = 1
        if padding:
            k = 2
        x = self.flux - m
        # complex FFT, since we need the conjugate
        s = fft.fft(x, k*n)
        ac = real(fft.ifft(s*s.conjugate()))
        ac /= ac[0]
        lag = self.time-self.time[0]
        if padding:
            return ac[:n], lag[:n]
        return ac[:n/2], lag[:n/2]
 
    def lowpass(self, low=10, decimate=False):
        """
        low pass Fourier filtering. 
        low = lower frequency limit (usually cycles per day)
        """
        from numpy import mean, fft, linspace, ones, searchsorted
        
        print 'lowpass on ',self.corotid
        raise(DeprecationWarning)
        m = mean(self.flux)
        
        n = len(self.flux)
        ft = fft.rfft(self.flux-m, n)
        freq = fft.fftfreq(n, d=self.dt())[:n/2]
        i = searchsorted(freq, low)
        i2 = i*2
        h = ones(n/2+1)
        h[i:i2] = linspace(1.0, 0.0, num=i2-i)
        h[i2:n] = 0.0
        ft *= h
        if decimate:
            n_org = n
            n = i2
            self.time = linspace(self.time[0], self.time[-1], num=n)
        
        # recalibrate amplitudes
        f = fft.irfft(ft, n)/(n_org/i2)
        self.flux = f+m 

    def resample(self, low=10):
        from scipy.signal import resample
        from numpy import array,linspace
        n = len(self.flux)
        num = int(n/((1./self.dt())/low))
        print n, num
        f = array([float(f) for f in self.flux])
        
        self.flux, self.time = resample(f, num, t=self.time)
        self.time = linspace(self.time[0], self.time[-1], num=num)

    def highpass(self, high=0.05):
        """
        low pass Fourier filtering. 
        low = lower frequency limit (usually cycles per day)
        """
        from numpy import mean, fft, zeros, linspace, searchsorted
        m = mean(self.flux)
        
        n = len(self.flux)
        n2 = 2*n
        ft = fft.rfft(self.flux-m, n2)
        freq = fft.fftfreq(n, d=self.dt())[:n-1]
        i = searchsorted(freq, high)
        i2 = i*2
        h = zeros(n+1)
        try:
            h[i:i2] = linspace(0.0, 1.0, num=i)
        except ValueError:
            print('h = %d, i = %d, i2 = %d' % (len(h),i,i2))
            exit()
        h[i2:n2] = 1.0
        ft *= h
        f = fft.irfft(ft)[:n] 
        #i = where(freq >= 0.0)
        self.flux = f+m 
    
    def butterworth(self):
        """
        low-pass filtering with butterworth but produces spikes! 
        """
        from scipy import signal
        from numpy import mean
        
        m = mean(self.flux)
        x = self.flux - m
        
        (N, Wn) = signal.buttord(wp=0.05, ws=0.5, gpass=2, gstop=30, analog=0)
        (b, a) = signal.butter(N, Wn, btype='low', analog=0, output='ba')
        
        # filtered output
        self.flux = signal.lfilter(b, a, x) + m                         

    def medfilter(self, length=11):
        """
        performs a median filter on the lightcurve
        """
        from scipy import signal
        self.flux = signal.medfilt(self.flux, length)
    
    def plot(self):
        """
        plot the lightcurve
        """
        import matplotlib.pyplot as plt
        plt.title = 'CorotID %s' % self.corotid
        plt.xtitle = 'HJD - 2455000'
        plt.ytitle = 'norm. flux'
        plt.grid()
        plt.plot(self.time-2455000.0, self.flux,'r')
        plt.show()
        
    def toasciifile(self, filename = 'lightcurve.txt'):
        f = open(filename,'w')
        for i in range(len(self.time)):
            f.write('%f\t%f\t%f\n' % (self.time[i], self.flux[i], self.mag[i]))
        f.close()
    
    def dt(self):
        """returns the time cadence of the lightcurve"""
        from numpy import median, diff
        
        return(median(diff(self.time)))  
    
    def jmpflt(self, sigma=1.4826):
        """
        removes jumps from the lightcurve
        performes sigma clipping on the absolute derivatives
        it is strongly advised to perform median filtering before applying 
        the filter and detrend after the filter. 
        """
        from numpy import diff, cumsum, array, std, sign
        
        def deriv(y, x):
            return (diff(y)/diff(x))
        def integrate(y, x):
            return (cumsum(y)*diff(x))

        x = self.time - round(self.time[0])
        y = self.flux

        y2 = deriv(y-1., x)
        stdlim = sigma*std(y2)
        
        y2 = [sign(i)*stdlim if (abs(i) > stdlim) else i for i in y2]
        
        ya = array(y2)
        self.time = self.time[:-1]
        self.flux = integrate(ya, x) + 1.
        
    def wavelets(self, noiseSigma = 16.0, wavefun=None, limit='soft'):
        """wavelet filtering of the lightcurve
        noiseSigma is in units of the std() of the lightcurve
        wavefun + one of [haar, db2, ...]
        limit is either, hard, soft, greater or less"""
        import pywt  # @UnresolvedImport
        from scipy import sqrt, log2
        from numpy import std, mean
        if wavefun is None:
            wavefun = 'sym5'
        wavelet = pywt.Wavelet(wavefun)

        n = len(self.flux)
        levels  = int( log2(n) ) 
        
        dec = pywt.wavedec(self.flux, wavelet, mode='cpd', level=levels)
        
        threshold = noiseSigma*std(self.flux)*sqrt(2*log2(n))
        
        if limit=='soft':
            nwc = map (lambda x: pywt.thresholding.soft(x, threshold), dec)
        elif limit=='hard':
            nwc = map (lambda x: pywt.thresholding.hard(x, threshold), dec)
        elif limit=='less': 
            nwc = map (lambda x: pywt.thresholding.less(x, threshold), dec)
        elif limit=='greater': 
            nwc = map (lambda x: pywt.thresholding.greater(x, threshold), dec)
                   
        filtered = pywt.waverec( nwc, wavelet, mode='cpd')
        
        if len(filtered) > len(self.flux):
            filtered = filtered[0:-1]
        dflux = self.flux - filtered
        self.flux = filtered
        self.flux = self.flux - mean(self.flux) + 1.
        return dflux
    
    def cut(self, fromtime, totime):
        """cut a certain time range from the lightcurve"""
        from numpy import nan, isnan, compress
        idx = self.time.searchsorted([fromtime,totime])
        self.time[idx[0]:idx[1]] = nan
        valid = ~isnan(self.time)
        self.time = compress(valid,self.time)
        self.flux = compress(valid,self.flux)

    def subtractscaled(self, times, flux):
        """subtract a scaled function of the given lightcurve in times and flux"""
        # interpolate to correct times
        from scipy.interpolate import interp1d
        f = interp1d(times, flux)
        if self.time[0]<times[0]:
            print str(self.time[0])+' < '+str(times[0])
        m = f(self.time)
        
        # build up matrix
        import numpy as np
        A = np.vstack([m, np.ones(len(m))]).T
        # get solution
        k,d = np.linalg.lstsq(A, self.flux)[0]

        #print k,d
        if k<0:
            print "Warning: k negative!"
        if d<0:
            print "Warning: d negative!"
        self.flux = self.flux - (k*m + d)
    
    def tofile(self, filename, overwrite=False):
        """save lightcurve to Fits-file"""
        import pyfits
        from numpy import zeros 
        c1 = pyfits.Column(name='DATEHEL', format='D', array=self.time-2451545.0)
        c2 = pyfits.Column(name='WHITEFLUX', format='D', array=self.flux)
        c3 = pyfits.Column(name='STATUS', format='J', array=zeros(len(self.flux)))
        tbhdu = pyfits.new_table([c1, c2, c3])
        hdu = pyfits.PrimaryHDU()
        
        hdu.header['COROTID'] = self.corotid
        hdu.header['RUN_CODE'] = self.hdr['RUN_CODE']
        hdu.header['HLFCCDID'] = self.hdr['HLFCCDID']
        hdu.header['WIN_ID'] = self.hdr['WIN_ID']
        hdu.header['MAGNIT_R'] = self.hdr['MAGNIT_R']
        hdu.header['MAGNIT_V'] = self.hdr['MAGNIT_V']
        hdu.header['MAGNIT_B'] = self.hdr['MAGNIT_B']
        
        thdulist = pyfits.HDUList([hdu, tbhdu])
        thdulist.writeto(filename, clobber=overwrite)
        