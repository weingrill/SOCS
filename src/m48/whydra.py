'''
Created on Aug 21, 2013

@author: jwe <jweingrill@aip.de>
'''
import config

class IsoChrone(dict):
    def __init__(self, filename = None):
        from numpy import loadtxt
        if filename is None:
            filename = config.datapath+'0p500Gyr_FeH0p0_Y0p277_AMLTsol.iso'
        a = loadtxt(filename)
        self['V'] = a[:,10]
        self['B-V'] = a[:,9]- a[:,10]

class WHydra(object):
    '''
    classdocs
    '''

    def __init__(self, field_name = 'M48field'):
        '''
        Constructor
        '''
        import ephem
        import datetime
        import pytz
        import astronomy as ast
        from datasource import DataSource

        self.field_name = field_name
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de') 
        
        kpno = ephem.Observer()
        #31.958036,-111.600578
        kpno.lon, kpno.lat = '-111.600578', '31.958036'
        kpno.horizon = '-0:34'
        kpno.elevation = 2096
        tzi = pytz.timezone('MST')
        #fmt = '%Y-%m-%d %H:%M:%S %Z%z'
        
        obsdate = datetime.datetime(2015,2,9,23,0,0, tzinfo=tzi)
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
        self.ebv = 0.031 # from Webda
        self.dm = 9.53 # from Webda

        target = {}
        self.table = []
        target['id'] = 6000
        target['name'] = 'M48center'
        target['ra'] = self.center[0]
        target['dec'] = self.center[1]
        target['class'] = 'C'
        self.table.append(target)
        self.targeted = []

        
    def _get_center(self):
        data = self.wifsip.query("""select avg(ra),avg(dec) from m48stars;""")
        print data[0][0],data[0][1]
        return (279.632565404658, 5.31001697525473)

    def priorities(self):
        """updates the priorieties in the m48stars table"""
        from scipy.interpolate import interp1d
        import numpy as np
        from functions import scaleto
        
        print 'calculate priorities ...'
        self.wifsip.execute("""update m48stars set priority=NULL;""")
        self.wifsip.commit()
        self.wifsip.execute("""update m48stars set priority=1.0 where vmag<16.5;""")
        self.wifsip.commit()
        
        data = self.wifsip.query("""SELECT tab, vmag 
                               FROM m48stars
                               WHERE not bv is null AND vmag<16.5
                               ORDER BY tab;""")
        tab = [d[0] for d in data]
        v = np.array([d[1] for d in data])
        p1 = scaleto(v,[1.0, 0.5])
        for i in range(len(tab)):
            print '%4d: %.3f --> %.3f' % (tab[i], v[i],p1[i])
            self.wifsip.execute("""UPDATE m48stars
                          SET priority = priority * %f
                          WHERE tab = %d;""" % (p1[i], tab[i]))
        self.wifsip.commit()   

        iso = IsoChrone('/work2/jwe/m48/data/output256520738433.dat')
        x = iso['V'] + self.dm
        y = iso['B-V'] + self.ebv
        i = np.argsort(x)
        print min(x),max(x)
        x = x[i]
        y = y[i]
        bvint = interp1d(x, y) #, kind='cubic'
        data = self.wifsip.query("""SELECT tab, vmag, bv 
                               FROM m48stars 
                               WHERE not bv is null AND V<16.5
                               ORDER BY tab;""")
        tab = [d[0] for d in data]
        v= np.array([d[1] for d in data])
        bv = np.array([d[2] for d in data])
        p = abs(bv - bvint(v))
        i = np.where(p > 0.7)
        p[i] = 0.7
        p = scaleto(p, [1.0, 0.0])
        for t in tab:
            i = tab.index(t)
            print '%d: V=%.3f B-V=%.3f c=%.3f p=%.3f' % (t,v[i],bv[i],bvint(v[i]),p[i])
            self.wifsip.execute("""UPDATE m48stars
                              SET priority = priority * %f
                              WHERE tab = %d;""" % (p[i], t))
        self.wifsip.commit()   

        data = self.wifsip.query("""SELECT tab, ra, dec
                               FROM m48stars
                               WHERE not ra is NULL AND not dec is NULL
                               ORDER BY TAB;""")
        tab = [d[0] for d in data]
        ra = np.array([d[1] for d in data])
        dec = np.array([d[2] for d in data])
        dist = np.sqrt((ra-self.center[0])**2+(dec-self.center[1])**2)
        i = np.where(dist > 0.5)
        #dist[i] = 0.5
        p = scaleto(dist, [1.0, 0.0])
        #p[i] = 0.0
        for t in tab:
            i = tab.index(t)
            print '%d: d=%.3f p=%.3f' % (t,dist[i],p[i])
            self.wifsip.execute("""UPDATE m48stars
                              SET priority = priority * %f
                              WHERE tab = %d;""" % (p[i], t))
        self.wifsip.commit()   
        self.wifsip.close()

    def setpointing(self, pointing):
        """set pointings according to hydrasim output!"""
         
        targets = ",".join([str(t) for t in self.targeted])
        print targets
        query = """UPDATE m48stars 
                   SET pointing=%d
                   WHERE tab in (%s)""" % (pointing,targets)
        self.wifsip.execute(query)
        self.wifsip.close()
        
    
    def from_database(self):
        
        # warum pointing NOT NULL?
        data = self.wifsip.query("""SELECT target, twomass, ra, dec
                               FROM wiyn2 
                               WHERE twomass IS NOT NULL
                               ORDER BY priority DESC;""")
        for d in data:
            target = {}
            target['id'] = int(d[0])
            target['name'] = str(d[1])
            target['ra'] = float(d[2])
            target['dec'] = float(d[3])
            target['class'] = 'O'
            self.table.append(target)

        fops = self.wifsip.query("""SELECT twomass,raj2000,dej2000 
                               FROM twomass 
                               WHERE twomass NOT IN 
                               (SELECT twomass 
                                FROM wiyn2 
                                WHERE twomass IS NOT NULL) 
                                AND circle'((%f, %f),0.5)' @> point(raj2000,dej2000) 
                                AND jmag>9.0 AND jmag<11.0
                               ORDER BY jmag 
                               LIMIT 2000;""" % self.center)
        for f in fops:
            target = {}
            target['id'] = 6001 + fops.index(f)
            target['name'] = str(f[0])
            target['ra'] = float(f[1])
            target['dec'] = float(f[2])
            target['class'] = 'F'
            self.table.append(target)

        #extra targets
        extra = self.wifsip.query("""SELECT target, twomass, ra, dec
                               FROM wiyn2 
                               WHERE pointing IS NULL
                               AND observed
                               AND circle'((%f, %f),0.5)' @> point(ra,dec) 
                               ORDER BY priority DESC;""" % self.center)
        for d in extra:
            target = {}
            target['id'] = int(d[0])
            target['name'] = str(d[1])
            target['ra'] = float(d[2])
            target['dec'] = float(d[3])
            target['class'] = 'E'
            self.table.append(target)

        self.wifsip.close()
    
    def skyfile(self, filename='/work1/jwe/Dropbox/IC4756/data/IC4756sky.coords'):
        import astronomy as ast
        f = open(filename,'r')
        lines = f.readlines()
        f.close()
        for l in lines:
            target = {}
            target['id'] = 9000 + lines.index(l)
            target['name'] = 'Sky'+str(target['id'])
            target['ra'] = ast.hms2dd(l[0:12])
            target['dec'] = ast.dms2dd(l[12:24])
            target['class'] = 'S'
            self.table.append(target)
    
    
    def tofile(self, filename= '/home/jwe/bin/hydra_simulator/whydra/IC4756field1.ast'):
        import astronomy as ast
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
            f.write(s)
        f.close()
        #s = 'pointing\tR.A.\tDec\tmag\tcomment'

    def fromfile(self, filename= '/home/jwe/bin/hydra_simulator/whydra/IC4756field1.hydra'):
        import astronomy as ast
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
                plt.scatter(t['ra'], t['dec'], marker='o', c='g')
            if t['class']=='E':
                plt.scatter(t['ra'], t['dec'], marker=',', c='gray')
            if t['class']=='F':
                plt.scatter(t['ra'], t['dec'], marker='^', c='r')
            if t['class']=='S':
                plt.scatter(t['ra'], t['dec'], marker='h', c='b')
            if t['class']=='C':
                plt.scatter(t['ra'], t['dec'], marker='+', c='k')
                
        plt.xlim([self.center[0]+0.55,self.center[0]-0.55])
        plt.ylim([self.center[1]-0.55,self.center[1]+0.55])
        plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.grid()
        #plt.show()
        plt.savefig(filename, format='pdf', dpi=300)
        plt.close()
    
    def dohydra(self, pointing):
        import subprocess
        subprocess.call(['/home/jwe/bin/hydra_simulator/dowhydra.sh',str(pointing)])     
        
"""
0000000000111111111122222222223333333333444444444455555555556666666666
0123456789012345678901234567890123456789012345678901234567890123456789    
9003 Sky9003              18 38 23.480 +05 31 55.20 S  STATUS=   8        
"""
if __name__ == '__main__':
    wh1 = WHydra('M48field1')
    wh1.priorities() 
    exit()
    wh1.from_database()
    wh1.skyfile()
    wh1.tofile()
    
    wh1.dohydra(1)
#     wh1.make_plot('/home/jwe/Downloads/IC4756field.pdf')
# 
    
# 
    path = '/home/jwe/bin/hydra_simulator/whydra/'
    for i in range(2,10):
        wh = WHydra('M48field%d' % i)
        wh.fromfile(path+'M48field%d.hydra' % (i-1))
        wh.setpointing(i-1)
        
        wh.tofile(path+'M48field%d.ast' % i)
        wh.dohydra(i)     
        del wh
  
