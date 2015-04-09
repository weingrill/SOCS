#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 9, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
    
class MilkyWay(object):
    '''
    classdocs
    '''


    def __init__(self, filename, magnitudes=False):
        '''
        Constructor
        '''
        try:
            self._load(filename)
        except IOError:
            print 'loading failed'
            self.query()
            self.make(magnitudes)
            self._save(filename)
        
    
    def query(self):
        from datasource import DataSource
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT radeg, dedeg, vt from tycho where vt<=15.193"""
        self.data = wifsip.query(query)
        print 'fetched %d stars' % len(self.data)
        
    def make(self, magnitudes=False):
        print 'filling matrix'
        self.z = np.zeros([180, 360])
        if magnitudes:
            for ra, dec, vt in self.data:
                x = np.trunc(ra)
                y = np.trunc(dec+90)
                # faintest star in tycho has 15.193 in vt
                z = 10.0**((15.193-vt)/2.5)
                try:
                    self.z[y,x] += z
                except ValueError:
                    print type(z),z
            self.z = -2.5*np.log10(self.z+1.0)
            self.z += 15.193
        else:
            for ra, dec in zip(self.x,self.y):
                x = np.trunc(ra)
                y = np.trunc((dec+90))
                self.z[y,x] += 1
    
    def _load(self, filename='/work2/jwe/tychomap.npy'):
        print 'loading matrix'
        self.z = np.load(filename)    
    
    def _save(self, filename='/work2/jwe/tychomap.npy'):
        print 'saving matrix'
        np.save(filename, self.z)

    def plot(self, show=True):
        """
        taken from http://matplotlib.org/1.4.3/examples/pylab_examples/contour_demo.html
        http://wiki.scipy.org/Cookbook/SignalSmooth
        """
        def gauss_kern(size, sizey=None):
            """ Returns a normalized 2D gauss kernel array for convolutions """
            from scipy import mgrid
            size = int(size)
            if not sizey:
                sizey = size
            else:
                sizey = int(sizey)
            x, y = mgrid[-size:size+1, -sizey:sizey+1]
            g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
            return g / g.sum()
        
        def blur_image(im, n, ny=None) :
            """ blurs the image by convolving with a gaussian kernel of typical
                size n. The optional keyword argument ny allows for a different
                size in the y direction.
            """
            from scipy import signal
            g = gauss_kern(n, sizey=ny)
            improc = signal.convolve(im,g, mode='same')
            return(improc)
        import matplotlib.pyplot as plt
        
        if show: plt.figure()
        #x = np.arange(0,360)
        x = np.arange(0, 24, 1./15)
        y = np.arange(-90, 90)
        print 'generating contour'
        z = blur_image(self.z, 9)
        print z.shape
        CS = plt.contour(x,y, z)
        if show: plt.show()
        return CS

if __name__ == '__main__':
    mw = MilkyWay('/work2/jwe/milkyway.npy', magnitudes=True)
    mw.plot()