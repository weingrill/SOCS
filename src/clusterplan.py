#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 19, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from datasource import DataSource
import ephem
from astropy.coordinates import SkyCoord  # @UnresolvedImport
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from tools import log
import datetime


class ClusterPlan(object):
    '''
    classdocs
    '''

    def __init__(self, clusterlist=None):
        '''
        Constructor:
        
        builds list of clusters either by criteria or by list given in clusterlist
        '''

        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')

        if clusterlist is None:
            query = self.selectclusters()
        else:
            if type(clusterlist) == str:
                clusterstring = "'%s'" % clusterlist
            else:
                clusterstring = ', '.join(["'%s'" % cluster for cluster in clusterlist])
            query = """SELECT name,ra,dec,diam,d,ebv,logage 
                FROM clusters
                WHERE name in (%s)""" % clusterstring
        print(clusterlist, query)

        result = self.wifsip.query(query)
        print('%d clusters found' % len(result))
        if len(result) == 0:
            raise ValueError('no clusters meet the criteria')

        self.current_year = datetime.datetime.now().year

        self.data = []
        for r in result:
            rec = {}
            rec['name'] = r[0]
            rec['ra'] = float(r[1])
            rec['dec'] = float(r[2])
            rec['diam'] = int(r[3])
            rec['d'] = int(r[4])
            rec['ebv'] = float(r[5])
            rec['age'] = round(10 ** float(r[6]) / 1e6, -1)  # convert to Myrs
            self.data.append(rec)

    def selectclusters(self):
        '''
        defines the criteria for the cluster selection and queries the database
        the result of the query is stored in the data property, which needs to be improved
        '''
        from math import log10
        criteria = {'minage': log10(125e6),
                    'maxage': log10(2600e6),
                    'maxebv': 0.31,
                    'maxd': 1500,
                    'mindec': -15.0,
                    'mindiam': 10,
                    'maxdiam': 80}

        query = """SELECT name,ra,dec,diam,d,ebv,logage from clusters 
            WHERE ((diam<=%(maxdiam)d and diam>=%(mindiam)d) or diam IS NULL) 
            AND logage>=%(minage)f and logage<=%(maxage)f
            AND dec>=%(mindec)f
            AND (name LIKE 'NGC%%' OR name LIKE 'IC%%')
            AND (ebv <= %(maxebv)f OR ebv IS NULL)
            AND (d <= %(maxd)d or d IS NULL)
            AND not observed
            order by ra""" % criteria

        return query

    def _eph2dt(self, ephemdate):
        """converts ephem.Date to datetime"""
        return datetime.datetime.strptime(str(ephemdate), "%Y/%m/%d %H:%M:%S")

    def plot(self):
        '''
        plots the milkyway density map and the selected clusters on the map.
        '''
        from milkyway import MilkyWay
        plt.figure(figsize=(10.69, 7.27))
        mw = MilkyWay('/work2/jwe/tychomap.npy', magnitudes=True)
        mw.plot(show=False)
        ra = [d['ra'] / 15. for d in self.data]
        dec = [d['dec'] for d in self.data]
        ebv = [d['ebv'] for d in self.data]
        name = [d['name'] for d in self.data]
        d = [d['diam'] * 4 for d in self.data]

        plt.scatter(ra, dec, edgecolor='none', c=ebv, s=d)
        plt.xlim(24, 0)
        plt.ylim(-15, 75)
        plt.minorticks_on()
        plt.grid()
        plt.xlabel('R.A.')
        plt.ylabel('Dec.')
        for r, d, n in zip(ra, dec, name):
            plt.text(r, d, n, fontsize=8)

        plt.savefig('/work2/jwe/SOCS/clusterplan.pdf', papertype='a4')
        plt.close()
        # plt.show()

    def list(self):
        '''
        lists the clusters that fulfill the criteria
        '''
        for d in self.data:
            print('%-15s %4dpc %4dMyr E(B-V)=%.2f %2d' % (d['name'], d['d'], d['age'], d['ebv'], d['diam']))
            # print self.time(d)

    def calc(self):
        """
        calculate the exposure for a solar like star
        """
        isofile = '/home/jwe/data/iso_01.000_Gyr.dat'

        a = np.loadtxt(isofile)
        iso_mass = a[:, 0]
        iso_mv = a[:, 5]
        iso_bv = a[:, 8]

        minbv, maxbv = 0.2, 1.45

        bv = np.arange(minbv, maxbv, 0.01)
        # mv = np.arange(min(iso_mv), max(iso_mv), 0.1)
        mass = np.arange(min(iso_mass), max(iso_mass), 0.01)

        i = np.where((iso_bv > minbv) & (iso_bv < maxbv) & (iso_mv > 2.0))
        p = np.polyfit(iso_bv[i], iso_mv[i], 5)
        # q = np.polyfit(iso_mv, iso_bv, 11)
        massmvp = np.polyfit(iso_mass, iso_mv, 12)
        massbvp = np.polyfit(iso_mass, iso_bv, 12)

        mv_int = np.polyval(massmvp, mass)
        bv_int = np.polyval(massbvp, mass)

        # p1 = [-1.33,  7.27,  0.75]
        print('[' + ', '.join(['%.2f' % pi for pi in p]) + ']')
        y = np.polyval(p, bv)
        print(np.polyval(p, [0.4, 1.4]))

        plt.scatter(iso_bv, iso_mv)

        plt.plot(bv, y, 'g')

        plt.plot(bv_int, mv_int, 'r')

        plt.axhline(4.83, color='y')
        plt.axvline(0.653, color='y')
        plt.xlabel('(B-V)')
        plt.ylabel('M$_V$')
        plt.ylim(plt.ylim()[::-1])
        plt.grid()
        plt.minorticks_on()
        plt.show()

    def loadfromfile(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
        dates = []
        hours = []
        darkhours = []

        for line in lines:
            split_line = line.split()
            date = ' '.join(split_line[0:2])
            dates.append(datetime.datetime.strptime(date, "%Y/%m/%d %H:%M:%S"))
            hours.append(float(split_line[2]))
            darkhours.append(float(split_line[3]))
        hours = np.array(hours)
        darkhours = np.array(darkhours)
        return dates, hours, darkhours

    def obstime(self):
        plt.figure(figsize=(29.6 / 2.54, 21. / 2.54))
        darkhours = np.zeros(365)
        for c in self.data:
            print(c['name'])
            filename = '/work2/jwe/SOCS/data/obstime %(name)s.txt' % c
            try:
                dates, hours, darkhours = self.loadfromfile(filename)
            except IOError:

                date0 = ephem.Date('%d/1/1 00:00:00' % self.current_year)
                hours = np.zeros(365)
                dates = []
                for day in range(365):
                    ephemdate = ephem.Date(date0 + day)
                    t = self.time(c, date=ephemdate)
                    dates.append(self._eph2dt(ephemdate))
                    hours[day] = t
                    if darkhours[day] == 0.0:
                        darkhours[day] = self.darktime(ephemdate)
                    log(filename, '%-20.20s %5.2f %5.2f %4.2f' % (ephemdate, t, darkhours[day], t / darkhours[day]))

            plt.plot(dates, hours, label=c['name'])
        plt.plot(dates, darkhours, 'k--')
        plt.grid()
        plt.yticks(np.arange(12))
        plt.ylabel('hours visible')
        plt.xlabel('date')
        xticks = [datetime.date(y, m, 1) for y in np.arange(self.current_year-1, self.current_year+1) for m in range(1, 13)]
        xlabels = [month.strftime('%b') for month in xticks]
        plt.xticks(xticks, xlabels)
        plt.xlim(dates[0], dates[-1])
        plt.title('clusterplan obstime')
        plt.legend(loc=9, fontsize='small')
        plt.savefig('/work2/jwe/SOCS/plots/clusterplan obstime %s.pdf' % c['name'])
        # plt.show()

    def darktime(self, date):
        izana = ephem.Observer()
        if date is None:
            date = ephem.Date('2015/8/2 00:00:00')
        izana.date = date
        izana.lat = '28.301195'
        izana.lon = '-16.509209'
        izana.horizon = '-0:34'
        # kpno.elevation = 2096
        sun = ephem.Sun(izana)  # @UndefinedVariable

        # set astronomical dawn
        izana.horizon = '-19:00'
        spset = izana.previous_setting(sun)
        snrise = izana.next_rising(sun)

        darktime = (ephem.Date(snrise) - ephem.Date(spset)) * 24.0
        return darktime

    def time(self, cluster, date=None, verbose=False):

        izana = ephem.Observer()
        if date is None:
            date = ephem.Date('2015/8/2 00:00:00')
        izana.date = date  # '2015/03/19 00:00:00'
        izana.lat = '28.301195'
        izana.lon = '-16.509209'
        izana.horizon = '-0:34'
        # izana.elevation = 2096
        sun = ephem.Sun(izana)  # @UndefinedVariable
        c = SkyCoord(cluster['ra'], cluster['dec'], frame='icrs', unit=(u.deg, u.deg))  # @UndefinedVariable
        rastr = c.ra.to_string(unit=u.hour, sep='::')  # @UndefinedVariable
        decstr = c.dec.to_string(unit=u.deg, sep='::')  # @UndefinedVariable
        ephemstr = '%s,f|O,%s,%s, 5.,2000' % (cluster['name'], rastr, decstr)
        clusterephem = ephem.readdb(ephemstr)
        eventlist = []
        # set airmass 2.0 
        izana.horizon = '30:00'
        clusterephem.compute(izana)
        cnrise = izana.next_rising(clusterephem)
        cnset = izana.next_setting(clusterephem)
        cprise = izana.previous_rising(clusterephem)
        cpset = izana.previous_setting(clusterephem)
        eventlist.append([self._eph2dt(cnrise), 'cluster rise (next)'])
        eventlist.append([self._eph2dt(cnset), 'cluster set (next)'])
        eventlist.append([self._eph2dt(cprise), 'cluster rise (prev)'])
        eventlist.append([self._eph2dt(cpset), 'cluster set (prev)'])

        # set astronomical dawn
        izana.horizon = '-19:00'
        spset = izana.previous_setting(sun)
        snrise = izana.next_rising(sun)
        eventlist.append([self._eph2dt(snrise), 'sunrise'])
        eventlist.append([self._eph2dt(spset), 'sunset'])

        eventlist.sort()
        sundown = False
        clusterup = False
        t0, t1 = None, None

        for t, e in eventlist:
            if verbose: print(t, e)
            if e == 'sunset':
                sundown = True
                if clusterup:
                    t0 = t

            if e == 'sunrise':
                sundown = False
                if clusterup:
                    t1 = t

            if e == 'cluster rise (next)' or \
                    e == 'cluster rise (prev)':
                clusterup = True
                if sundown: t0 = t
            if e == 'cluster set (next)' or \
                    e == 'cluster set (prev)':
                clusterup = False
                if sundown: t1 = t

        if verbose: print(t0, t1)
        if t0 is None and t1 is None:
            return 0.0
        return (ephem.Date(t1) - ephem.Date(t0)) * 24.0


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='SOCS cluster observation planning')
    parser.add_argument('--plot', action='store_true', help='plot clusters')
    parser.add_argument('-obstime', action='store_true', help='plot observation schedule')
    parser.add_argument('-list', action='store_true', help='list clusters')
    parser.add_argument('-calc', action='store_true', help='calculate exposure time')
    parser.add_argument('clusters', help='clusters to calculate')

    args = parser.parse_args()

    cp = ClusterPlan(clusterlist=args.clusters)
    if args.plot: cp.plot()
    if args.obstime: cp.obstime()
    if args.list: cp.list()
    if args.calc: cp.calc()
