#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 17, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datasource import DataSource
# from datetime import timedelta
from astronomy import jd
import datetime


class ObservedClusters(object):
    '''
    produces a timeline of observations for a given object
    '''

    def __init__(self):
        '''
        Constructor
        '''
        wifsip = DataSource(database='stella',
                            user='stella',
                            host='pera.aip.de')

        query = """SELECT distinct object 
        FROM frames 
        WHERE object LIKE 'NGC%rot NW' OR object LIKE 'M%rot NW' OR object LIKE 'NGC%rot'
        """

        result = wifsip.query(query)
        self.objects = [r[0].rstrip(' NW') for r in result]

        self.data = []
        for obj in self.objects:
            query = """select distinct floor(jd) "startdate", expt "duration" 
            from frames 
            where object like '%s%%'
            order by startdate
            """ % obj
            result = wifsip.query(query)

            startdates = np.array([r[0] for r in result])
            # durations = np.array([timedelta(seconds=r[1]) for r in result])
            durations = np.array([r[1] for r in result])
            rec = {'object': obj, 'startdates': startdates, 'durations': durations}
            self.data.append(rec)

        # self.startdates = np.array([r[1] for r in result])
        # self.enddates = np.array([r[2] for r in result])
        print('%d records' % len(result))

    def plot(self, show=False):

        plt.figure(figsize=(25 / 2.54, 18 / 2.54))
        # plt.plot_date(x=self.dates, y=self.count)
        # plt.title('%s (total: %d frames)' % (self.target, sum(self.count)))
        plt.title('Observed open clusters produced by SOCS/observed.py')
        plt.grid(which='both')
        # plt.ylabel('number of frames')
        # plt.minorticks_on()

        for data in self.data:
            obj = data['object']
            print(obj)
            n = len(data['startdates'])
            bottom = np.array([self.objects.index(obj)] * n)
            left = data['startdates']

            width = np.ones(n) * 0.5

            plt.barh(bottom, width, left=left, linewidth=0)
        bottom = np.arange(len(self.objects))
        clusters = [obj.rstrip(' rot') for obj in self.objects]
        plt.yticks(bottom + 0.4, clusters)

        dates = [datetime.datetime(y, m, 1) for y in range(2013, 2018) for m in range(1, 13)]
        datenames = [date.strftime('%b %Y') for date in dates]
        xticks = [jd(date) for date in dates]
        i0 = datenames.index('Sep 2013')
        i1 = datenames.index('Dec 2017')
        plt.xticks(xticks[i0:i1], datenames[i0:i1], rotation=45)
        plt.xlim(xticks[i0], xticks[i1])
        if show:
            plt.show()
        else:
            plt.savefig('/home/jwe/Downloads/observed clusters.pdf')
        plt.close()


if __name__ == '__main__':
    os = ObservedClusters()
    os.plot(show=False)
