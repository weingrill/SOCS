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
            WHERE diam<40 and diam>15 
            AND logage>8.477 and logage<9.3
            AND dec>-6
            AND (name like 'NGC%' or name like 'IC%')
            AND ebv < 0.5
            AND d < 1500"""
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

    def plot(self):
        import matplotlib.pyplot as plt
        ra = [d['ra']/15. for d in self.data]
        dec = [d['dec'] for d in self.data]
        ebv = [d['ebv'] for d in self.data]
        name = [d['name'] for d in self.data]
        d = [d['diam']*3 for d in self.data]
        plt.scatter(ra, dec, edgecolor='none', c=ebv, s=d)
        plt.xlim(24,0)
        for r,d,n in zip(ra, dec, name):
            plt.text(r, d, n)
        plt.show()
        
    def list(self):
        for d in self.data[-4:-1]:
            print '%-15s %4dpc %4dMyr %.2f %2d' % (d['name'],d['d'],d['age'],d['ebv'],d['diam'])
            self.time(d)
        
    def time(self, cluster):
        import ephem
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        
        izana = ephem.Observer()
        date = ephem.Date('2015/03/19 00:00:00')
        izana.date = '2015/03/19 00:00:00'
        izana.lat = '28.301195'
        izana.lon = '-16.509209'
        izana.horizon = '-0:34'
        #izana.elevation = 2096
        sun = ephem.Sun()  # @UndefinedVariable
        
        c = SkyCoord(cluster['ra']*15, cluster['dec'], 'icrs', unit=(u.deg, u.deg))  # @UndefinedVariable
        rastr = str(c.ra).replace('d', ':').replace('m', ':').rstrip('s')
        decstr = str(c.dec).replace('d', ':').replace('m', ':').rstrip('s')
        ephemstr = '%s,f|O,%s,%s, 5.,2000' % (cluster['name'],rastr,decstr)
        
        clusterephem = ephem.readdb(ephemstr)
        
        # set airmass 2.0 
        izana.horizon = '30:00'
        clusterephem.compute(izana)
        cnrise = izana.next_rising(clusterephem)
        cnset =  izana.next_setting(clusterephem)
        cprise = izana.previous_rising(clusterephem)
        cpset =  izana.previous_setting(clusterephem)
        #print cnrise
        print 'objset ', cprise
        print 'objrise', cnset
        #print cpset
        #set astronomical dawn
        izana.horizon = '-15:00'
        spset =  izana.previous_setting(sun)
        snset =  izana.next_setting(sun)
        snrise = izana.next_rising(sun)
        sprise = izana.previous_rising(sun)
        #print sprise
        #print snset
        print 'sunset ', spset
        print 'sunrise', snrise
        if cprise<spset: cstart=spset 
        else: cstart=cprise
        if cnset>snrise: cstop=snrise
        else: cstop=cnset  
        print (cstop-cstart)*24
        return
        for i in range(3):
            newdate = date + i
            izana.date = newdate
            sun.compute(izana)
            at = izana.next_antitransit(sun)
            
            izana.date = at 
            print izana.sidereal_time(), izana.date 
        

if __name__ == '__main__':
    cp = ClusterPlan()
    #cp.plot()
    cp.list()
    
        
    
