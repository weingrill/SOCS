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
        
        self.target = target
        print self.target
        wifsip = DataSource(database='wifsip', 
                           user='sro', 
                           host='pina.aip.de')
        
        query = """SELECT date_trunc('day', datesend- interval '12 hour') "datum", count(objid)
        FROM frames 
        WHERE object like '%s%%'
        AND 1.*matched/stars>0.8
        AND backgrnd<500
        AND moondist>40.0
        GROUP BY datum
        ORDER BY datum;""" % (self.target)
        
        result = wifsip.query(query)
        
        self.dates = [r[0] for r in result]
        self.count = [r[1] for r in result]
        
    def plot(self, show=False):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from numpy import sum
        fig = plt.figure(figsize=(10,6))
        #plt.plot_date(x=self.dates, y=self.count)
        plt.bar(self.dates,self.count)
        fig.autofmt_xdate()
        plt.title('%s (total: %d frames)' % (self.target, sum(self.count)))
        plt.grid(which='both')
        plt.ylabel('number of frames')
        plt.minorticks_on()
        
        if show:
            plt.show()
        else:
            plt.savefig('/home/jwe/Downloads/%s timeline.pdf' % self.target)
        plt.close()
    
if __name__ == '__main__':
    import sys
    if len(sys.argv)>1:
        for cluster in sys.argv[1:]:
            tl = TimeLine(cluster)
            tl.plot(False)
    else:
        for cluster in ('NGC 6940 rot','NGC 2281 rot', 'M 48 rot', 'NGC 6633 rot',
                        'NGC 1528 rot', 'M 67 rot'):
            tl = TimeLine(cluster)
            tl.plot(False)
