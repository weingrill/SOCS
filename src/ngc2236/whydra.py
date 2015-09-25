#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 21, 2013

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
hydrasimpath = '/home/jwe/bin/hydra_simulator/'
hydrapath = '/home/jwe/bin/hydra_simulator/whydra/'

class IsoChrone(dict):
    def __init__(self, filename = None):
        from numpy import loadtxt
        
        a = loadtxt(filename)
        self['V'] = a[:,10]
        self['B-V'] = a[:,9]- a[:,10]

class WHydra(object):
    '''
    classdocs
    '''

    def __init__(self, field_name = 'NGC2236field'):
        '''
        Constructor
        '''
        import ephem
        import datetime
        import pytz
        import astronomy as ast
        from datasource import DataSource
        from cluster import Cluster
        
        self.field_name = field_name
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost) 
        self.corot = DataSource(database='corot', user='sro', host='oldpina.aip.de') 
        
        kpno = ephem.Observer()
        #31.958036,-111.600578
        kpno.lon, kpno.lat = '-111.600578', '31.958036'
        kpno.horizon = '-0:34'
        kpno.elevation = 2096
        tzi = pytz.timezone('MST')
        #fmt = '%Y-%m-%d %H:%M:%S %Z%z'
        
        obsdate = datetime.datetime(2016,2,9,23,0,0, tzinfo=tzi)
        kpno.date = obsdate + datetime.timedelta(7*ephem.hour)
        d = kpno.date.datetime()
        
           
        self.input_epoch = 2000
        self.current_epoch = d.year + d.timetuple().tm_yday/365.
        
        # LST at mid-exposure in decimal hours
        lst = kpno.sidereal_time()
        print 'local siderial time: ',lst
        self.siderial_time = ast.hms2hh(str(lst))
        self.exposure_length = 40.0/60.0 # in hours
        self.wavelength = 6000 # in Angstroem
        self.cable = 'BLUE'
        self.weighting = 'STRONG'
        self.guidewave = 6000
        
        self.center = self._get_center()
        
        c = Cluster('NGC 2236')
        self.ebv = 0.55 #c['ebv'] # ?= 0.28
        self.dm = 13.8 #ast.distance_modulus(c['d'])

        target = {}
        self.table = []
        target['id'] = 6000
        target['name'] = 'NGC2236center'
        target['ra'] = self.center[0]
        target['dec'] = self.center[1]
        target['class'] = 'C'
        self.table.append(target)
        self.targeted = []
        # number of sky fibers
        self.skies = 6
        # number of field orientation probes
        self.fops = 6
        self.tweakcenter = False
        self.tweakangle = False
        
        print 'center', self.center
        print 'E(B-V)', self.ebv
        print 'DM', self.dm
        
    def _get_center(self):
        data = self.wifsip.query("""select avg(ra),avg(dec) from ngc2236;""")
        return (data[0][0],data[0][1])

    def priorities(self, verbose = False):
        """updates the priorities in the ngc2236 table"""
        import numpy as np
        from functions import scaleto
        
        def makeplot(bv, v, p, filename=None, isobv=None, isov=None):
            """plot the priorities"""
            from matplotlib import pyplot
            pyplot.scatter(bv, v, c=p, edgecolor='none', alpha=0.75)
            if not isobv is None:
                pyplot.plot(isobv, isov, 'k')
            pyplot.xlim(-0.2,2.0)
            pyplot.ylim(16.5,7.5)
            pyplot.title('NGC 2236')
            pyplot.xlabel('B - V')
            pyplot.ylabel('V mag')
            
            if filename is None: 
                pyplot.show()    
            else:
                pyplot.savefig(filename, dpi=300)
            pyplot.close()

        def plotradec(ra, dec, p, filename=None):
            """plot the priorities"""
            from matplotlib import pyplot
            pyplot.scatter(ra, dec, c=p, edgecolor='none', alpha=0.75)
            pyplot.title('NGC 2236')
            pyplot.xlabel('R.A.')
            pyplot.ylabel('Dec')
            
            if filename is None: 
                pyplot.show()    
            else:
                pyplot.savefig(filename, dpi=300)
            pyplot.close()

            
        print 'calculate priorities ...'
        self.wifsip.execute("UPDATE ngc2236 SET priority=NULL, pointing=NULL;")
        self.wifsip.execute("""UPDATE ngc2236 
            SET priority=1.0 
            WHERE vmag<16.5 
            AND NOT tab IS NULL;""")
        
        data = self.wifsip.query("""SELECT tab, vmag, bv 
                               FROM ngc2236
                               WHERE not bv is null AND vmag<16.5
                               ORDER BY tab;""")
        tab = [d[0] for d in data]
        v = np.array([d[1] for d in data])
        bv = np.array([d[2] for d in data])
        p1 = scaleto(v,[1.0, 0.8])
        makeplot(bv, v, p1, filename=config.plotpath+'vmag_priorities.pdf')
        
        print len(tab),'stars brighter V<16.5'
        for i in range(len(tab)):
            if verbose: print '%4d: %.3f --> %.3f' % (tab[i], v[i],p1[i])
            self.wifsip.execute("""UPDATE ngc2236
                          SET priority = priority * %f
                          WHERE tab = %d;""" % (p1[i], tab[i]), commit=False)
        self.wifsip.commit()   

        lengood = self.wifsip.query("select count(starid) from ngc2236 where good;")
        print lengood[0][0],'stars with periods'
        
        self.wifsip.execute("""UPDATE ngc2236
                          SET priority = priority * 0.8
                          WHERE NOT good;""")
        self.wifsip.commit()   

        self.wifsip.execute("""UPDATE ngc2236
                          SET priority = priority * 0.5
                          WHERE period IS NULL;""")
        self.wifsip.commit()   

        iso = IsoChrone(config.datapath+'output885516794937.dat')
        x = iso['V'] + self.dm
        y = iso['B-V'] + self.ebv
        data = self.wifsip.query("""SELECT tab, vmag, bv 
                               FROM ngc2236 
                               WHERE not bv is null AND vmag<16.5
                               ORDER BY tab;""")
        tab = [d[0] for d in data]
        v= np.array([d[1] for d in data])
        bv = np.array([d[2] for d in data])
        p = np.zeros(len(tab))
        print len(tab),'stars for isochrone priority'
        
        for i in range(len(tab)):
            p[i] = np.min(abs(x-v[i])+abs(y-bv[i]))
            
        #p[p > 0.4] = 0.4
        p = scaleto(p, [1.0, 0.1])
        makeplot(bv, v, p, filename=config.plotpath+'iso_priorities.pdf', isobv=y, isov=x)
        for t in tab:
            i = tab.index(t)
            if verbose: 
                print '%d: V=%.3f B-V=%.3f p=%.3f' % (t, v[i], bv[i], p[i])
            self.wifsip.execute("""UPDATE ngc2236
                              SET priority = priority * %f
                              WHERE tab = %d;""" % (p[i], t), commit=False)
        self.wifsip.commit()   
        
        data = self.wifsip.query("""SELECT tab, ra, dec
                               FROM ngc2236
                               WHERE not ra is NULL AND not dec is NULL
                               AND priority > 0.0
                               ORDER BY TAB;""")
        tab = [d[0] for d in data]
        print len(tab),'stars for distance priority'
        ra = np.array([d[1] for d in data])
        dec = np.array([d[2] for d in data])
        dist = np.sqrt((ra-self.center[0])**2+(dec-self.center[1])**2)
        from functions import gauss
        p = gauss(dist, 1.0, 0.0, 7.0/60.)
        #p = scaleto(dist, [1.0, 0.0])
        plotradec(ra, dec, p, filename=config.plotpath+'coord_priorities.pdf')
        for t in tab:
            i = tab.index(t)
            try:
                if verbose:
                    print '%d: d=%.3f p=%.3f' % (t,dist[i],p[i])
                self.wifsip.execute("""UPDATE ngc2236
                              SET priority = priority * %f
                              WHERE tab = %d;""" % (p[i], t), commit=False)
            except TypeError:
                print t
        self.wifsip.commit() 
        
        data = self.wifsip.query("""SELECT bv,vmag, priority
                               FROM ngc2236
                               WHERE priority > 0.0
                               ORDER BY tab;""")
        
        bv = np.array([d[0] for d in data])
        vmag = np.array([d[1] for d in data])
        p = np.array([d[2] for d in data])
        makeplot(bv,vmag, p, filename=config.plotpath+'priorities.pdf')
          
    def setpointing(self, pointing):
        """set pointings according to hydrasim output"""
        #TODO: to be redone! 
        targets = ",".join([str(t) for t in self.targeted])
        print targets
        query = """UPDATE ngc2236 
                   SET pointing=%d
                   WHERE tab in (%s)""" % (pointing,targets)
        self.wifsip.execute(query)
    
    def from_database(self, maglimit=16.5):
        """
        load targets from database that have not priority set yet and are
        within the maglimit
        """
        data = self.wifsip.query("""SELECT tab, starid, ra, dec, vmag
                               FROM ngc2236 
                               WHERE priority > 0.0
                               AND vmag<%f
                               AND pointing IS NULL
                               ORDER BY priority DESC;""" % maglimit)
        for d in data:
            target = {}
            target['id'] = int(d[0])
            target['name'] = 'Star%04dm%.2f' % (int(d[0]),float(d[4]))
            target['ra'] = float(d[2])
            target['dec'] = float(d[3])
            target['class'] = 'O'
            self.table.append(target)
        """
        hydrawiynmanual.pdf p.41:
        The stars selected for use by the FOPs should fall in the magnitude 
        range 10<V<14. If possible, keep the range in magnitude of your FOPs 
        sample in each field to 3 magnitudes or less so that the intensity of 
        each star falls within the dynamic range of the FOPs TV camera. 
        Including the FOPs star magnitudes in the configuration file may
        also be useful later when setting up the field at the telescope.
        """
        query = """SELECT -2.1*log(lc_meang)+21.984 "vmag", alpha, delta, corotid 
        FROM corot 
        WHERE run_code='SRa02' AND hlfccdid like 'E1%%' 
        AND point(alpha,delta) <@ circle(point(%f, %f),0.5) 
        AND NOT lc_meang IS NULL AND lc_meang>19000 
        ORDER BY lc_meang DESC;""" % self.center
        fops = self.corot.query(query)
        print len(fops),'FOPs stars'
        for f in fops:
            target = {}
            i = fops.index(f)
            target['id'] = 6001 + i
            target['name'] = 'FOP%dm%.2f' % (int(f[3]),float(f[0]))
            target['ra'] = float(f[1])
            target['dec'] = float(f[2])
            target['class'] = 'F'
            self.table.append(target)

        #extra targets
        params = dict(zip(['ra','dec'],self.center))
        query="""select corotid, alpha,delta
        from corot where run_code='SRa02' 
        and hlfccdid like 'E1%%'
        and point(alpha,delta) <@ circle(point(%(ra)f, %(dec)f),0.5)
        order by point(alpha,delta) <-> point(%(ra)f, %(dec)f)
        LIMIT 999;"""% params
        
        extra = self.corot.query(query)
        print len(extra),'extra stars'
        for d in extra:
            target = {}
            target['id'] = 7000+extra.index(d)
            target['name'] = str(d[0])
            target['ra'] = float(d[1])
            target['dec'] = float(d[2])
            target['class'] = 'E'
            self.table.append(target)
    
    def skyfile(self, filename='/work2/jwe/NGC2236/skyfile.txt'):
        """
        load the sky positions
        """
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        f = open(filename,'r')
        lines = f.readlines()
        f.close()
        print len(lines),'sky positions'
        
        for l in lines:
            try:
                c = SkyCoord(l, 'icrs', unit=(u.hourangle, u.deg))  # @UndefinedVariable
            except ValueError:
                print 'empty line'
            else:
                target = {}
                target['id'] = 9000 + lines.index(l)
                target['name'] = 'Sky'+str(target['id'])
                target['ra'] = c.ra.degree
                target['dec'] = c.dec.degree
                target['class'] = 'S'
                self.table.append(target)
    
    def tofile(self, filename = None, verbose=False):
        """
        writes the .ast file that finally goes into whydra for processing
        """
        import astronomy as ast
        if filename is None:
            filename = hydrapath+self.field_name+'.ast'
        f = open(filename,'w')
        f.write('FIELD NAME: %s\n' % self.field_name)
        f.write('INPUT EPOCH: %.2f\n' % self.input_epoch)
        f.write('CURRENT EPOCH: %.2f\n' % self.current_epoch)
        f.write('SIDERIAL TIME: %.3f\n' % self.siderial_time)
        f.write('EXPOSURE LENGTH: %.2f\n' % self.exposure_length)
        f.write('WAVELENGTH: %.0f.\n' % self.wavelength)
        f.write('CABLE: %s\n' % self.cable)
        f.write('WEIGHTING: %s\n' % self.weighting)
        f.write('GUIDEWAVE: %.0f.\n' % self.guidewave)
        f.write('#0000000011111111112222222222333333333344444444445555\n')                          
        f.write('#2345678901234567890123456789012345678901234567890123\n')                          
        for t in self.table:
            ra = ast.dd2hms(t['ra'])
            dec = ast.dd2dms(t['dec'])
            s = '%04d %-20s %02d %02d %06.3f %+02.2d %02d %05.2f %1s\n' % \
                (t['id'],t['name'],ra[0],ra[1],ra[2], dec[0], dec[1], dec[2], t['class'])
            if verbose: print s.rstrip('\n')
            f.write(s)
        f.close()
        #s = 'pointing\tR.A.\tDec\tmag\tcomment'

    def fromfile(self, filename = None):
        import astronomy as ast
        if filename is None:
            filename = hydrapath+self.field_name+'.hydra'
        
        f = open(filename,'r')
        lines = f.readlines()
        f.close()
        for l in lines:
            if l[0]=='#':
                l = ''
            #if l.find('FIELD NAME:')>=0:
            # self.field_name = l[l.find(':')+1:].strip()
               
            if l.find('INPUT EPOCH:')>=0:
                self.input_epoch = float(l[l.find(':')+1:])
            if l.find('CURRENT EPOCH:')>=0:
                self.current_epoch = float(l[l.find(':')+1:])
            if l.find('SIDERIAL TIME:')>=0:
                self.siderial_time = float(l[l.find(':')+1:])
            if l.find('EXPOSURE LENGTH:')>=0:
                self.exposure_length = float(l[l.find(':')+1:])
            if l.find('WAVELENGTH:')>=0:
                self.wavelength = float(l[l.find(':')+1:])
            if l.find('CABLE:')>=0:
                self.cable = l[l.find(':')+1:].strip()
            if l.find('WEIGHTING:')>=0:
                self.weighting = l[l.find(':')+1:].strip()
            if l.find('GUIDEWAVE:')>=0:
                self.guidewave = float(l[l.find(':')+1:])
            if len(l.strip())>=53:
                target = {
                'id': int(l[0:4]),
                'name': l[5:25],
                'ra': ast.hms2dd(l[26:38]),
                'dec': ast.dms2dd(l[39:51]),
                'class': l[52]}
                if l.find('STATUS=')>=55:
                    status=l[55:66].strip()
                    if target['class'] in ['O','E']:
                        if (status=='STATUS=OK' or status=='STATUS=EDGE'):
                            self.table.append(target)
                        else:
                            self.targeted.append(target['id'])
                    elif target['class'] in ['C','F','S']:
                        self.table.append(target)
                        
    def make_plot(self, filename):
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=[18/2.54,18/2.54])
        fov = plt.Circle((self.center[0], self.center[1]), 0.5, facecolor='none', edgecolor='g')
        fig = plt.gcf()
        fig.gca().add_artist(fov)
#        plt.scatter(ra,dec, marker=',', c='k', s=0.1)
        
        
        for t in self.table:
            if t['class']=='O':
                plt.scatter(t['ra'], t['dec'], marker='o', c='g', edgecolor='none')
            if t['class']=='E':
                plt.scatter(t['ra'], t['dec'], marker=',', c='gray', edgecolor='none')
            if t['class']=='F':
                plt.scatter(t['ra'], t['dec'], marker='^', c='r', edgecolor='none')
            if t['class']=='S':
                plt.scatter(t['ra'], t['dec'], marker='h', c='b', edgecolor='none')
            if t['class']=='C':
                plt.scatter(t['ra'], t['dec'], marker='+', c='k', edgecolor='none')
                
        plt.xlim([self.center[0]+0.55,self.center[0]-0.55])
        plt.ylim([self.center[1]-0.55,self.center[1]+0.55])
        plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.grid()
        #plt.show()
        plt.savefig(filename, dpi=300)
        plt.close()

        #data = self.wifsip.query("""SELECT vmag, bv, priority 
        #                       FROM m48stars 
        #                       WHERE priority>0.1;""")
        #v= np.array([d[0] for d in data])
        #bv = np.array([d[1] for d in data])
        #p = [d[2] for d in data]
        #from matplotlib import pyplot
        #pyplot.scatter(bv, v, c=p)
        #pyplot.xlim(-0.1,1.6)
        #pyplot.ylim(16.5,7.5)
        #pyplot.show()    
        #plt.savefig(filename, dpi=300)
        #plt.close()

    def dohydra(self):
        """
        write the commands file and
        execute the shell script
        """
        import subprocess
        
        f = open('/home/jwe/bin/hydra_simulator/cmds.%s' % self.field_name,'wt')
        f.write('%s.ast\n' % self.field_name)
        f.write('%s\n' % self.field_name)
        f.write('%d\n' % self.fops)
        f.write('%d\n' % self.skies)
        if self.tweakangle:
            f.write('y\n')
        else:
            f.write('n\n')
        if self.tweakcenter:
            f.write('y\n')
        else:
            f.write('n\n')
        f.close()
        
        subprocess.call(['/home/jwe/bin/hydra_simulator/dowhydra.sh',self.field_name])  
        
    def getconcentricities(self):
        """
        fetches the current concentricities file from the WIYN web page
        """
        import urllib2
        response = urllib2.urlopen('http://www.wiyn.org/concentricities')
        html = response.read()
        f = open(hydrapath+'concentricities','wt')
        f.write(html)
        f.close()
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NGC2236 WHydra preparation')
    parser.add_argument('-p', '--priorities', action='store_true', 
                        help='calculate priorities')
    parser.add_argument('-c', '--concentricities', action='store_true', 
                        help='get concentricities file')
    parser.add_argument('--plot', action='store_true', 
                        help='make the plot')
    parser.add_argument('--run', action='store_true', 
                        help='calculate the whole run')

    parser.add_argument('--bright', action='store_true', 
                        help='make bright pointing')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='verbose output')

    wh0 = WHydra('NGC2236field0')
    args = parser.parse_args()
    if args.priorities: wh0.priorities() 
    if args.concentricities: wh0.getconcentricities()
    if args.bright:
        wh0.from_database(maglimit=12.0)
        wh0.skyfile()
        wh0.tofile()
        wh0.dohydra()
        wh0.fromfile()
        wh0.setpointing(0)
    if args.run:
        wh1 = WHydra('NGC2236field1')
        wh1.from_database()
        wh1.skyfile()
        wh1.tofile()
        wh1.dohydra()
        del wh1
        for i in range(2,10):
            wh = WHydra('NGC2236field%d' % i)
            wh.fromfile(hydrapath+'NGC2236field%d.hydra' % (i-1))
            wh.setpointing(i-1)
            
            wh.tofile()
            wh.dohydra()     
            del wh

    if args.plot: wh0.make_plot('/home/jwe/Downloads/NGC2236field.pdf')
      
