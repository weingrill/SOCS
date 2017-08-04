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
        if target == 'M 48':
            datelimit = 'TRUE'
        elif target == 'NGC 6940':
            datelimit = "datesend > '2017-01-01'"
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
        
        todates = self.dates[1:] 
        fromdates = self.dates[0:-1]
        datedeltas = [todate - fromdate for todate, fromdate in zip(todates,fromdates)]
        
        # find gaps larger than 10 days
        offsets = [0]
        for i,datedelta in enumerate(datedeltas):
            if datedelta> datetime.timedelta(days=10):
                offsets.append(i)
                offsets.append(i+1)
                
        if i not in offsets:
            offsets.append(i)        
        
        print offsets
        
        if len(offsets) == 2:
            self.length = self.dates[-1] - self.dates[0]
            start = 0
            end = len(self.dates)
        else:
            self.length = self.dates[1] - self.dates[0]
            for o1, o2 in zip(offsets[1:], offsets[:-1]):
                length = self.dates[o1] - self.dates[o2] 
                print length.days, self.dates[o1], self.dates[o2], np.sum(self.count[o2:o1]) 
                if length.days > self.length.days and np.sum(self.count[o2:o1])>10:
                    self.length = length
                    start = o2
                    end = o1
        print start,end            
        self.dates = self.dates[start:end]
        self.count = self.count[start:end]    
        self.cumsum = np.cumsum(self.count)
        
        
    def plot(self, show=False):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        import datetime
        
        fig = plt.figure(figsize=(10,6))
        plt.plot(self.dates, self.cumsum, drawstyle='steps-mid')
        fig.autofmt_xdate()
        plt.title('%s (usable: %d frames, %d days)' % (self.target, self.cumsum[-1], self.length.days))
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
