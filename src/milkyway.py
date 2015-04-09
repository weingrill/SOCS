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


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT radeg, dedeg, vt from tycho WHERE vt<=15.193"""
        result = self.wifsip.query(query)
        self.x = np.array([r[0] for r in result])
        self.y = np.array([r[1] for r in result])
        self.z = np.zeros([180,360])
        for ra, dec, vt in result:
            x = np.trunc(ra)
            y = np.trunc(dec+90)
            # faintest star in tycho has 15.193 in vt
            z = 10.0**((15.193-vt)/2.5)
            self.z[y,x] += z
        self.z = -2.5*np.log10(self.z)
        
        np.save('/work2/jwe/tychomap', self.z)
        
    def plot(self):
        """
        taken from http://matplotlib.org/1.4.3/examples/pylab_examples/contour_demo.html
        """
        import matplotlib
        #import numpy as np
        import matplotlib.cm as cm
        import matplotlib.mlab as mlab
        import matplotlib.pyplot as plt
        
        plt.figure()
        CS = plt.contour(self.z)
        plt.show()

if __name__ == '__main__':
    mw = MilkyWay()
    mw.plot()