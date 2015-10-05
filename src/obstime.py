'''
Created on 11.02.2015

@author: Joerg Weingrill <jweingrill@gmail.com>
'''

def clusterephem(clustername):
    import cluster
    c = cluster.Cluster('NGC 2236') 
    #print c.coordinatestring
    coords = ','.join(c.coordinatestring)
    coords = coords.replace(' ', ':')
    mag = c.solarmagnitude
    diam = c['diam']*60
    result = '%s,f|O,%s,%.1f,2000,%i' % (clustername, coords, mag,diam)
    return result 

if __name__ == '__main__':
    import ephem, datetime
    from math import pi
    
    kpno = ephem.Observer()
    kpno.date = '2016/02/13 21:00:00' # UT!
    kpno.lon, kpno.lat = '-111.600578', '31.958036'
    #31.958035, -111.600585
    kpno.horizon = '20'
    kpno.elevation = 2096

    delta =  datetime.timedelta(7*ephem.hour) 
    #ic4756 = ephem.readdb('IC 4756,f|O,18:39: 0,+05:27, 5.,2000,3120')
    print clusterephem('NGC 2236')
    clu = ephem.readdb(clusterephem('NGC 2236'))
    
    clu.compute(kpno)
    rise = ephem.localtime(kpno.next_rising(clu, use_center=True))
    transit = ephem.localtime(kpno.next_transit(clu))
    clu_set = ephem.localtime(kpno.next_setting(clu))

    sun = ephem.Sun()  # @UndefinedVariable
    moon = ephem.Moon()
    moon.compute(kpno)

    sep =  ephem.separation((clu.az, clu.alt), (moon.az, moon.alt))

    #We relocate the horizon to get twilight times
    kpno.horizon = '-18' #-6=civil twilight, -12=nautical, -18=astronomical
    sun.compute(kpno)
    beg_twilight=ephem.localtime(kpno.next_setting(sun, use_center=True)) #Begin civil twilight
    end_twilight=ephem.localtime(kpno.next_rising(sun, use_center=True, start=beg_twilight)) #End civil twilight

    print moon
    print '%.1f, %.1f  - %.1f, %.1f' % (clu.az * 180.0 / pi, clu.alt * 180.0 / pi, moon.az * 180.0 / pi, moon.alt * 180.0 / pi)
    print 'sep: %.1f deg' % (sep * 180.0 / pi)
    print 'dark     ', (beg_twilight-delta).strftime("%Y-%m-%d %H:%M")
    print 'rise:    ', (rise-delta).strftime("%Y-%m-%d %H:%M")
    print 'transit: ', (transit-delta).strftime("%Y-%m-%d %H:%M")
    print 'set:     ', (clu_set-delta).strftime("%Y-%m-%d %H:%M")
    print 'end dark ', (end_twilight-delta).strftime("%Y-%m-%d %H:%M")
    if rise>beg_twilight:
        begin = rise
    else:
        begin = beg_twilight
    if clu_set<end_twilight:
        end = clu_set
    else:
        end = end_twilight 
    print begin-delta,'-',end-delta
    print 'obstime', end-begin
