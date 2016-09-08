#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 14, 2014

@author: jwe
'''

class TimeLine(object):
    '''
    produces a timeline of observations for a given object
    '''


    def __init__(self, target):
        '''
        Constructor
        '''
        from datasource import DataSource
        import datetime
        import numpy as np
        
        self.target = target
        print self.target
        wifsip = DataSource(database='stella', 
                           user='stella', 
                           host='pera.aip.de')
        #SELECT date_trunc('day', datesend- interval '12 hour') "datum", count(objid)
        if target == 'M 67':
            datelimit = "datesend < '2015-07-01'"
        elif target == 'NGC 2281':
            datelimit = "datesend > '2016-06-01'"
        elif target == 'M 48':
            datelimit = "datesend < '2014-06-01'"
        else:
            datelimit = 'TRUE'
        
        params = {'target': self.target, 'datelimit': datelimit}
        query = """SELECT datesend "datum", 1
         FROM frames  
         WHERE (object like '%(target)s rot NW' or object like '%(target)s rot') 
         AND expt>60 AND filter='V'
         AND %(datelimit)s
         ORDER BY datum;""" % params
        
        result = wifsip.query(query)
        
        self.dates = [r[0] for r in result]
        self.count = [r[1] for r in result]
        
        """base = self.dates[0]
        datedelta = self.dates[-1] - self.dates[0]
        numdays = datedelta.days+1
        date_list = [base + datetime.timedelta(days=x) for x in range(0, numdays)]
        counts = [0]*numdays
        for i,d in enumerate(self.dates):
            j = date_list.index(d)
            c = self.count[i]
            counts[j] = c
            print i, d, j, c
        
        self.dates = date_list
        self.count = counts
        
        for d, c in zip(self.dates, self.count):
            print d, c
        """
        
        self.count = np.cumsum(self.count)    
        
    def plot(self, show=False):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        import datetime
        
        fig = plt.figure(figsize=(10,6))
        plt.plot(self.dates, self.count, drawstyle='steps-mid')
        fig.autofmt_xdate()
        plt.title('%s (total: %d frames)' % (self.target, self.count[-1]))
        plt.grid(which='both')
        plt.ylabel('cumulative frames per night')
        xticks = [datetime.date(y, m, 1) for y in np.arange(2012, 2020) for m in range(1,13)]
        xlabels = [month.strftime('%B %Y') for month in xticks]
        plt.xticks(xticks, xlabels)
        
        plt.xlim(self.dates[0],self.dates[-1])
        
        if show:
            plt.show()
        else:
            plt.savefig('/work2/jwe/SOCS/plots/%s timeline.pdf' % self.target)
        plt.close()
    
if __name__ == '__main__':
    import sys
    if len(sys.argv)>1:
        for cluster in sys.argv[1:]:
            tl = TimeLine(cluster)
            tl.plot(False)
    else:
        for cluster in ('NGC 2422', 
                        'NGC 2323', 
                        'M 48', 
                        'NGC 6709',
                        'NGC 2281', 
                        'M 67', 
                        'NGC 1528', 
                        'NGC 6940',
                        'NGC 6633'):
            tl = TimeLine(cluster)
            tl.plot(False)
