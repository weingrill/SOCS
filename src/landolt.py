#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 6, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np


def medfit(x, y):
    from numpy import sum
    from scipy.optimize import minimize
    def fun(p):
        return sum(abs(y - p[0] - p[1] * x))

    p0 = [0.0, 1.0]

    res = minimize(fun, p0, method='Nelder-Mead')
    # print res
    return res.x


class Landolt(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.landolt = DataSource(database='stella',
                                  host='pera',
                                  user='stella')

    def match(self, objid):

        queryparam = {'objid': objid}
        query = """SELECT name, vmag, mag_auto, bv, magerr_auto, flux_auto
            FROM frames, phot, landolt
            WHERE phot.objid = '%(objid)s'
            AND filter='V'
            AND frames.objid=phot.objid
            AND circle(phot.coord,0.4/3600.) @> circle(landolt.coord,0)
            AND phot.flags<4;""" % queryparam
        result = self.landolt.query(query)
        self.names = [r[0] for r in result]
        self.vmag = np.array([r[1] for r in result])
        self.mag_auto = np.array([r[2] for r in result])
        self.bv = np.array([r[3] for r in result])
        self.err = np.array([r[4] for r in result])
        self.flux = np.array([r[5] for r in result])
        import matplotlib.pyplot as plt
        x = self.bv
        y = self.vmag - self.mag_auto
        n = self.names
        p = medfit(x, y)
        plt.title('corr = %.3f; k" = %.3f' % (p[0], p[1]))
        # plt.plot(self.bv, self.vmag-self.mag_auto, 'o')
        plt.errorbar(x, y, self.err, fmt='go')
        x1 = np.array(plt.xlim())
        plt.plot(x1, p[0] + p[1] * x1, 'r--')
        for xi, yi, ni in zip(x, y, n):
            plt.text(xi, yi, ni)
        plt.show()
        # print result

    def extinction(self):
        import pickle
        try:
            with open('/work1/jwe/Landolt/data/20130328A.pickle', 'rb') as picklefile:
                result = pickle.load(picklefile)
        except:
            query = """SELECT name, airmass, mag_auto, flux_auto
            FROM frames, phot, landolt
            WHERE filter = 'V'
            AND frames.objid like '20130328A%%'
            AND frames.objid = phot.objid
            AND circle(phot.coord,0.4/3600.) @> circle(landolt.coord,0)
            AND phot.flags<4
            ORDER BY name;"""
            result = self.landolt.query(query)
            for r in result:
                print(r)
            with open('/work1/jwe/Landolt/data/20130328A.pickle', 'wb') as picklefile:
                pickle.dump(result, picklefile)

        self.names = [r[0] for r in result]
        self.airmass = np.array([r[1] for r in result])
        self.mag_auto = np.array([r[2] for r in result])
        self.flux = np.array([r[3] for r in result])
        self.mag = -2.5 * np.log10(self.flux)

        mags = []
        airmasses = []
        n = len(result)
        for i in range(n):
            for j in range(i + 1, n):
                if self.airmass[i] > self.airmass[j]:
                    mags.append(self.mag[i] - self.mag[j])
                    airmasses.append(self.airmass[i] - self.airmass[j])
        mags = np.array(mags)
        airmasses = np.array(airmasses)

        import matplotlib.pyplot as plt

        k = np.median(mags[mags > -.1] / airmasses[mags > -.1])
        print(k)

        plt.scatter(airmasses, mags, edgecolor='none', alpha=0.75)
        # plt.xlim(0.0, 4.0)
        plt.ylim(-1.0, 1.0)
        x1 = np.array(plt.xlim())
        plt.plot(x1, k * x1, 'r-')
        plt.xlabel('$\Delta$X')
        plt.show()


if __name__ == '__main__':
    l = Landolt()
    l.extinction()
    # l.match('20130327A-0077-0003')
