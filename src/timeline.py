#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 14, 2014

@author: jwe
'''
from datasource import DataSource
import datetime
import numpy as np

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

gap_length =10


class TimeLine(object):
    '''
    produces a timeline of observations for a given object
    '''

    def __init__(self, target, frot=False):
        '''
        Constructor
        '''
        self.target = target

        wifsip = DataSource(database='stella',
                            user='stella',
                            host='pera.aip.de')
        if target == 'M 48':
            date_limit = 'TRUE'
        elif target == 'NGC 6940':
            date_limit = "datesend > '2017-01-01'"
        else:
            date_limit = 'TRUE'

        if frot:
            expt = 6
            self.survey = 'frot'
        else:
            expt = 60
            self.survey = 'rot'

        params = {'target': self.target,
                  'date_limit': date_limit,
                  'expt': expt,
                  "survey": self.survey,
                  "frot": frot}
        query = 'SELECT datesend "datum", 1 ' \
                "FROM frames " \
                "WHERE (object like '%(target)s %(survey)s NW' or object like '%(target)s %(survey)s')" \
                " AND expt>%(expt)d AND filter='V'" \
                " AND %(date_limit)s " \
                "ORDER BY datum;" % params

        result = wifsip.query(query)

        if len(result)==0:
            print("No observations found")
            return
        self.dates = [r[0] for r in result]
        self.count = [r[1] for r in result]

        todates = self.dates[1:]
        fromdates = self.dates[0:-1]
        datedeltas = [todate - fromdate for todate, fromdate in zip(todates, fromdates)]

        # find gaps larger than 10 days
        offsets = [0]
        i = 0
        for i, datedelta in enumerate(datedeltas):
            if datedelta > datetime.timedelta(days=gap_length):
                offsets.append(i)
                offsets.append(i + 1)

        if i not in offsets:
            offsets.append(i)

        start = 0
        end = len(self.dates)

        if len(offsets) == 2:
            self.length = self.dates[-1] - self.dates[0]
            start = 0
            end = len(self.dates)
        else:
            self.length = self.dates[1] - self.dates[0]
            for o1, o2 in zip(offsets[1:], offsets[:-1]):
                length = self.dates[o1] - self.dates[o2]
                print(length.days, self.dates[o1], self.dates[o2], np.sum(self.count[o2:o1]))
                if length.days > self.length.days and np.sum(self.count[o2:o1]) > gap_length:
                    self.length = length
                    start = o2
                    end = o1
        print(start, end)
        self.dates = self.dates[start:end]
        self.count = self.count[start:end]
        self.cumsum = np.cumsum(self.count)

    def plot(self, show=False):
        fig = plt.figure(figsize=(10, 6))
        plt.plot(self.dates, self.cumsum, drawstyle='steps-mid')
        fig.autofmt_xdate()
        plt.title('%s %s (usable: %d frames, %d days)' % (self.target, self.survey, self.cumsum[-1], self.length.days))
        plt.grid(which='both')
        plt.ylabel('cumulative frames per night')
        xticks = [datetime.date(y, m, 1) for y in np.arange(2012, 2020) for m in range(1, 13)]
        xlabels = [month.strftime('1. %B %Y') for month in xticks]
        plt.xticks(xticks, xlabels)

        plt.xlim(self.dates[0], self.dates[-1])

        if show:
            plt.show()
        else:
            plt.savefig('/work2/jwe/SOCS/plots/%s %s timeline.pdf' % (self.target, self.survey))
        plt.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='SOCS cluster timeline calculation')
    parser.add_argument('-frot', action='store_true', help='plot frot timeline')
    parser.add_argument('-show', action='store_true', default=False, help='plot frot timeline')
    parser.add_argument('cluster', help='cluster to calculate')

    args = parser.parse_args()
    tl = TimeLine(args.cluster, frot=args.frot)
    tl.plot(args.show)
