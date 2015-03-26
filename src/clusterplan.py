#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Mar 19, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class ClusterPlan(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT name,ra,dec,diam,d,ebv,logage from clusters 
            WHERE diam<60 and diam>15 
            AND logage>8.477 and logage<9.46
            AND dec>-6
            AND (name like 'NGC%' or name like 'IC%')
            AND ebv < 0.5
            AND d < 1500
            AND not observed"""
        result = self.wifsip.query(query)
        self.data = []
        for r in result:
            rec = {}
            rec['name'] = r[0]
            rec['ra'] = float(r[1])
            rec['dec'] = float(r[2])
            rec['diam'] = int(r[3])
            rec['d'] = int(r[4])
            rec['ebv'] = float(r[5])
            rec['age'] = round(10**float(r[6])/1e6,-1)
            self.data.append(rec)
    
    def _eph2dt(self, ephemdate):
        """converts ephem.Date to datetime"""
        import datetime
        return datetime.datetime.strptime( str(ephemdate), "%Y/%m/%d %H:%M:%S" )

    def plot(self):
        import matplotlib.pyplot as plt
        ra = [d['ra']/15. for d in self.data]
        dec = [d['dec'] for d in self.data]
        ebv = [d['ebv'] for d in self.data]
        name = [d['name'] for d in self.data]
        d = [d['diam']*4 for d in self.data]
        plt.scatter(ra, dec, edgecolor='none', c=ebv, s=d)
        plt.xlim(24,0)
        plt.xlabel('R.A.')
        plt.ylabel('Dec.')
        for r,d,n in zip(ra, dec, name):
            plt.text(r, d, n)
        
        plt.savefig('/work2/jwe/SOCS/clusterplan.pdf')
        plt.close()
        #plt.show()
        
    def list(self):
        for d in self.data[-3:]:
            print '%-15s %4dpc %4dMyr %.2f %2d' % (d['name'],d['d'],d['age'],d['ebv'],d['diam'])
            print self.time(d)
    
    def obstime(self):
        import ephem
        import numpy as np
        import matplotlib.pyplot as plt
        
        darkhours = np.zeros(365)
        for c in self.data:
            print c['name']
            date0 = ephem.Date('2015/1/1 00:00:00')
            hours = np.zeros(365)
            dates = []
            for day in range(365):
                ephemdate = ephem.Date(date0 + day)
                t = self.time(c, date = ephemdate)
                print ephemdate, t
                dates.append(self._eph2dt(ephemdate))
                hours[day] = t
                if darkhours[day] == 0.0:
                    darkhours[day] = self.darktime(ephemdate) 
                
            
            plt.plot(dates,hours, label=c['name'])
        plt.plot(dates, darkhours, 'k--')
        plt.grid()
        plt.minorticks_on()
        plt.legend(loc=9, fontsize='small')
        plt.show()
    
    def darktime(self, date):
        import ephem
        
        izana = ephem.Observer()
        if date is None:
            date = ephem.Date('2015/8/2 00:00:00')
        izana.date = date #'2015/03/19 00:00:00'
        izana.lat = '28.301195'
        izana.lon = '-16.509209'
        izana.horizon = '-0:34'
        #izana.elevation = 2096
        sun = ephem.Sun(izana)  # @UndefinedVariable
        
        #set astronomical dawn
        izana.horizon = '-19:00'
        spset =  izana.previous_setting(sun)
        snrise = izana.next_rising(sun)
        
        darktime = (ephem.Date(snrise)-ephem.Date(spset))*24.0
        return darktime
        
    def time(self, cluster, date = None, verbose = False):
        import ephem
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        
        
        izana = ephem.Observer()
        if date is None:
            date = ephem.Date('2015/8/2 00:00:00')
        izana.date = date #'2015/03/19 00:00:00'
        izana.lat = '28.301195'
        izana.lon = '-16.509209'
        izana.horizon = '-0:34'
        #izana.elevation = 2096
        sun = ephem.Sun(izana)  # @UndefinedVariable
        c = SkyCoord(cluster['ra'], cluster['dec'], 'icrs', unit=(u.deg, u.deg))  # @UndefinedVariable
        rastr = c.ra.to_string(unit=u.hour, sep='::') # @UndefinedVariable
        decstr = c.dec.to_string(unit=u.deg,sep='::') # @UndefinedVariable
        ephemstr = '%s,f|O,%s,%s, 5.,2000' % (cluster['name'],rastr,decstr)
        clusterephem = ephem.readdb(ephemstr)
        eventlist = []
        # set airmass 2.0 
        izana.horizon = '30:00'
        clusterephem.compute(izana)
        cnrise = izana.next_rising(clusterephem)
        cnset =  izana.next_setting(clusterephem)
        cprise = izana.previous_rising(clusterephem)
        cpset =  izana.previous_setting(clusterephem)
        eventlist.append([self._eph2dt(cnrise),'cluster rise (next)'])
        eventlist.append([self._eph2dt(cnset),'cluster set (next)'])
        eventlist.append([self._eph2dt(cprise),'cluster rise (prev)'])
        eventlist.append([self._eph2dt(cpset),'cluster set (prev)'])
        
        #set astronomical dawn
        izana.horizon = '-19:00'
        spset =  izana.previous_setting(sun)
        snrise = izana.next_rising(sun)
        eventlist.append([self._eph2dt(snrise),'sunrise'])
        eventlist.append([self._eph2dt(spset),'sunset'])
        
        
        eventlist.sort()
        sundown = False
        clusterup = False
        t0, t1 = None, None
        
        for t,e in eventlist:
            if verbose: print t,e
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
                
        if verbose: print t0,t1
        if t0 is None and t1 is None:
            return 0.0
        return (ephem.Date(t1)-ephem.Date(t0))*24.0

if __name__ == '__main__':
    cp = ClusterPlan()
    #cp.plot()
    cp.obstime()
    #cp.list()
    
        
    
