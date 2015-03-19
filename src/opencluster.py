'''
Created on Oct 3, 2013

@author: jwe <jweingrill@aip.de>
'''
class OpenCluster(object):
    '''
    classdocs
    '''

    def __init__(self, objectname='', uname=None, ra=None, dec=None, obsmode=None):
        """
        initializes the OpenCluster Object with certain values
        obsmode can be 'rot' for rotation photometry
                    or 'cmd' for BVI photometry
        """
        
        self.objectname = objectname
        self.startdate = '2013-03-21'
        self.enddate = '2018-12-31'
        self.priority = 1.0
        self.telescope = 'STELLA-I'
        self.withfocus = True
        self.withacquire = True
        self.withguiding = True
        self.guiderchoice = 'piggy-back'
        self.title = 'SOCS'
        if uname is None:
            self.uname = '%s %s' % (objectname, obsmode)
        else:
            self.uname = uname
        self.propid = 'cluster.survey'
        self.abstract = 'Photometric monitoring of open stellar clusters'
        self.pi = 'Weingrill'
        self.affil = 'CORE'
        self.team = 'Barnes'
        self.mode = {'mode': 'Clusters'}
        self.camera = {'camera':'direct',
                       'XOffCCD':0,
                       'YOffCCD':0,
                       'XSizeCCD':2050,
                       'YSizeCCD':4100,
                       'XBinCCD':1,
                       'YBinCCD':1}
        self.sequence = {'sequence':'FullFilters',
                        'offset':0.0}
        self.object = {'ObjectName':objectname,
                       'RA':ra,
                       'Dec':dec}
        self.constraints = {'MoonDistance.Min':40,
                            'SolHeight.Max':-16.0,
                            'AirmassTarget.Max':2.0,
                            'AltTarget.Min':30.0}
        self.file = ''
        self.filename = ''
        self.fields = 1
        self.obsmode = obsmode
        
        if not (obsmode in ['BVI','BVR','uvby','Hby','rot']):
            raise TypeError('%s not a valid observation mode' % obsmode)

        if obsmode == 'BVI':
            self.sequence['ExposureTime'] = 20.0 
            self.sequence['ExposureRepeat'] = 6
            self.sequence['ExposureIncrease'] = '1,2,4,10,15,30'
            self.sequence['FilterSequence'] = 'I,V,B,I,V,B'
            self.mode['mode'] = 'Clusters'
            self.mode['pickdelay'] = 0
            self.mode['pernight'] = 2 
            self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
            # zerofraction is the length of one exposure in days
            self.mode['zerofraction'] = 2.0/24.0
            self.mode['impact'] = 1.0

        if obsmode == 'BVR':
            self.sequence['ExposureTime'] = 20.0 
            self.sequence['ExposureRepeat'] = 6
            self.sequence['ExposureIncrease'] = '1,2,4,10,15,30'
            self.sequence['FilterSequence'] = 'R,V,B,R,V,B'
            self.mode['mode'] = 'Clusters'
            self.mode['pernight'] = 2 
            self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
            # zerofraction is the length of one exposure in days
            self.mode['zerofraction'] = 2.0/24.0
            self.mode['impact'] = 1.0

        if obsmode == 'uvby':
            self.sequence['ExposureTime'] = 30.0 
            self.sequence['ExposureRepeat'] = 8
            self.sequence['ExposureIncrease'] = '2,3,4,5,8,12,16,20'
            self.sequence['FilterSequence'] = 'y,b,v,u,y,b,v,u'
            self.mode['mode'] = 'Clusters'
            self.mode['pernight'] = 2 
            self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
            #self.mode['timeout'] = self.timeout
            # zerofraction is the length of one exposure in days
            self.mode['zerofraction'] = 2.0/24.0
            self.mode['impact'] = 1.0

        if obsmode == 'Hby':
            self.sequence['ExposureTime'] = 60.0 
            self.sequence['ExposureRepeat'] = 6
            self.sequence['ExposureIncrease'] = '1,1.5,1,3,2,6'
            self.sequence['FilterSequence'] = 'y,b,hbw,hbn,haw,han'
            self.mode['mode'] = 'Clusters'
            self.mode['pernight'] = 2 
            self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
            #self.mode['timeout'] = self.timeout
            # zerofraction is the length of one exposure in days
            self.mode['zerofraction'] = 2.0/24.0
            self.mode['impact'] = 1.0

        if obsmode == 'rot':
            self.sequence['ExposureTime'] = 30.0
            self.sequence['ExposureRepeat'] = 3
            self.sequence['ExposureIncrease'] = '2,10,20'
            self.sequence['FilterSequence'] = 'V,V,R'
            self.mode['mode'] = 'Clusters'
            #self.mode['timeout'] = self.timeout # duration*fields*1000
            self.mode['pernight'] = 4 # can be refined
            self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
            # zerofraction is the length of one exposure in days
            self.mode['zerofraction'] = 1.0/24.0
            self.mode['impact'] = 1.0
        
        if self.object['RA'] is None or self.object['Dec'] is None:
            self.get_coordiantes()
    
    def plot_ephem(self, obsdate=None):
        import matplotlib.pyplot as plt
        
        import ephem
        import datetime
        from numpy import pi,empty
        
        from astronomy import airmass
        
        stella = ephem.Observer()
        #stella.lon, stella.lat = '13.104659', '52.404963' # Potsdam
        stella.lat, stella.lon = 31.9583, -111.59867 # KPNO
        stella.lat, stella.lon = 28.301214,-16.509246
        sun, moon = ephem.Sun(), ephem.Moon()  # @UndefinedVariable
        
        stella.pressure = 0
        stella.horizon = '-0:34'
        stella.elevation = 2000
        
        ephemstr = ','.join([self.objectname,
                             'f|O',
                             self.ra_str,
                             self.dec_str,
                             '5.0'])
        
        ocluster = ephem.readdb(ephemstr)
        ocluster.compute()
        
        print 'Moonrise:', stella.previous_rising(moon)
        print 'Moonset: ', stella.next_setting(moon)
        print 'Sunrise: ', stella.previous_rising(sun)
        print 'Sunset:  ', stella.next_setting(sun)
        
        if obsdate is None:
            today = datetime.datetime.today()
        else: today = datetime.datetime.strptime(obsdate,'%Y/%m/%d %H:%M:%S')
            
        #dt =  datetime.timedelta(days=14)
        #today += dt
        
        sun_alt = empty(24)
        moon_alt = empty(24)
        hours = range(24)
        ocluster_alt = empty(24)
        for h in hours:
            today = today.replace(hour=h,minute=0,second=0)
            stella.date = today 
            sun.compute(stella)
            moon.compute(stella)
            ocluster.compute(stella)
            sun_alt[h] = float(sun.alt)
            moon_alt[h] = float(moon.alt)
            ocluster_alt[h] = float(ocluster.alt)
            #print alt[h]
        
        fig = plt.figure()
        ax_h = fig.add_subplot(111)
        
        ax_h.set_ylim(0,90) 
        ax_h.set_xlim(0,24) 
    
        ax_airmass = ax_h.twinx()
        
        ax_h.set_xticks(hours)
        heights = ax_h.get_yticks()
        am = airmass(heights)
        aml = ['%.2f ' % a for a in am]
        ax_airmass.set_ylim(0.,90.)
        ax_airmass.set_yticklabels(aml)
        ax_h.grid()
        ax_h.plot(hours, sun_alt*180.0/pi,'yo')
        ax_h.plot(hours, moon_alt*180.0/pi,'go')
        ax_h.plot(hours, ocluster_alt*180.0/pi,'k')
        
        ax_h.set_xlabel("hours")
        ax_h.set_ylabel("height (degrees)")
        ax_airmass.set_ylabel("airmass")
        
        plt.draw()  
        plt.show()            
        
    def get_coordiantes(self):
        """
        queries the object coordinates
        """
        import PySimbad
        import astronomy as ast
        
        self.coords = PySimbad.SimbadCoord(self.object['ObjectName'])
        
        if self.coords.find('+')>0:
            p = self.coords.find('+')
        else:
            p = self.coords.find('-')
        ra = self.coords[:p].strip()    
        dec = self.coords[p:].strip()
        self.ra_str = ra
        self.dec_str = dec
        self.object['RA'] = ast.hms2dd(ra)
        self.object['Dec'] = ast.dms2dd(dec)
        self.coords=[self.object['RA'],self.object['Dec']]
        
    def plan_wifsip(self, nfields=4):
        """
        returns new subframes
        """
        self.fields = nfields
        d = 1320.2/3600.0 # 1320.2arcsec is the fov of WiFSIP
        d2 = d/2
        cra, cdec = self.object['RA'],self.object['Dec']
        if nfields == 4:
            names = ['NW','NE','SW','SE']
            fields = [(cra - d2, cdec + d2),
                 (cra + d2, cdec + d2),
                 (cra - d2, cdec - d2),
                 (cra + d2, cdec - d2)]
        if nfields == 5:
            names = ['C','NW','NE','SW','SE']
            fields = [(cra    , cdec),
                 (cra - d2, cdec + d2),
                 (cra + d2, cdec + d2),
                 (cra - d2, cdec - d2),
                 (cra + d2, cdec - d2)]
        
        subframes = []
        for f in fields:
            i = fields.index(f)
            subframes.append(OpenCluster(objectname=self.object['ObjectName'],
                                         uname = '%s %s' % (self.uname, names[i]), 
                                         ra = f[0], dec = f[1], obsmode = self.obsmode))
        #do some deep copy of the object:
        for sf in subframes:
            sf.fields = self.fields
            sf.startdate = self.startdate
            sf.enddate = self.enddate
            sf.telescope = self.telescope
            sf.withfocus = self.withfocus
            sf.withacquire = self.withacquire
            sf.withguiding = self.withguiding
            sf.title = self.title
            #sf.uname = self.uname # is set at init
            sf.propid = self.propid
            sf.abstract = self.abstract
            sf.pi = self.pi
            sf.affil = self.affil
            sf.team = self.team
            sf.title = self.title
            sf.mode = self.mode
            sf.camera = self.camera
            sf.sequence = self.sequence
            #sf.object = self.object # new coordinates and object mustn't change!
            sf.constraints = self.constraints
            
        return subframes
            
    def plot(self, axis=None):
        """
        plot the fov
        """
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(1.)
        self.tycho()
        ra, dec = self.object['RA'],self.object['Dec']
        d2 = 0.5*1320.2/3600.0
        ras = [ra-d2, ra+d2, ra+d2, ra-d2, ra-d2]
        das = [dec-d2, dec-d2, dec+d2, dec+d2, dec-d2]
        if axis is None:
            plt.plot(ras, das)
        else:
            axis.plot(ras, das)
        plt.show()
    
    @property
    def exposuretime(self):
        return self.sequence['ExposureTime']
    
    @property
    def pernight(self):
        return self.mode['pernight']

    @property
    def exposurerepeat(self):
        seq = self.sequence['ExposureIncrease'] #'5*1,5*1.5,5*3'
        sseq = seq.split(',')
        repeat = 0
        for s in sseq:
            if s.find('*')>0:
                repeat += int(s.split('*')[0])
            else:
                repeat += 1
        return repeat

    @property
    def duration(self):
        """
        getter for duration
        """
        expt = self.exposuretime
        seq = self.sequence['ExposureIncrease'] #'5*1,5*1.5,5*3'
        sseq = seq.split(',')
        repeat = 0
        duration = 0.0
        for s in sseq:
            if s.find('*')>0:
                repeat += int(s.split('*')[0])
                duration += expt*int(s.split('*')[0])*float(s.split('*')[1])
            else:
                repeat += 1
                duration += expt*float(s)
        if self.sequence['ExposureRepeat']<>repeat:
            print "Warning: ExposureRepeat different!"      
        self.sequence['ExposureRepeat'] = repeat
        return duration
    
    @property
    def timeout(self):
        """calculate timeout"""
        return 0.0
        #return self.duration * self.fields * 1000.0
    
    @property
    def zerofraction(self):
        return self.duration/86400.0
    
    @property
    def pickdelay(self):
        return self.duration * self.fields
    
    def tofile(self, path = './'):
        """
        stores the data to a save file
        """
        #ngc6940_sw_bvi_20130926.save
        import os
        from datetime import date, datetime
        self.file = self.uname.lower().replace(' ','_')+'.xml'
        filename, _ =  os.path.splitext(self.file)
        todaystr = date.strftime(date.today(),'%Y%m%d')
        self.filename = os.path.join(path,filename + '_' + todaystr + '.save')
        def str2bin(boolean):
            '''return the binary value to an appropriate string'''
            return str(boolean).lower()
        
        f = open(self.filename, 'wt')
        f.write('#User input\n')
        #Thu Sep 26 10:33:48 CEST 2013
        f.write('#%s\n' % datetime.strftime(datetime.now(),'%a %b %d %H:%M:%S %Z %Y'))
        f.write('startdate=%s\n' % self.startdate)
        f.write('enddate=%s\n' % self.enddate)
        f.write('priority=%.1f\n' % self.priority)
        f.write('TELESCOPE=%s\n' % self.telescope)
        f.write('withfocus=%s\n' % str2bin(self.withfocus))
        f.write('withacquire=%s\n' % str2bin(self.withacquire))
        f.write('withguiding=%s\n' % str2bin(self.withguiding))
        f.write('title=%s\n' % str(self.title))
        f.write('uname=%s\n' % str(self.uname))
        f.write('propid=%s\n' % self.propid)
        f.write('abstract=%s\n' % self.abstract)
        f.write('pi=%s\n' % self.pi)
        f.write('affil=%s\n' % self.affil)
        f.write('team=%s\n' % self.team)
        f.write('mode=%s\n' % self.mode['mode'])

        if self.mode['mode']=='Clusters':
            f.write('mode.pickdelay=%d\n' % self.pickdelay)
            f.write('mode.pernight=%d\n' % self.mode['pernight'])
            f.write('mode.period_day=%f\n' % self.mode['period_day'])
            f.write('mode.zerofraction=%f\n' % self.zerofraction)
            f.write('mode.impact=%f\n' % self.mode['impact'])
            
        f.write('camera=%s\n' % self.camera['camera'])
        f.write('camera.XOffCCD=%d\n' % self.camera['XOffCCD'])
        f.write('camera.YOffCCD=%d\n' % self.camera['YOffCCD'])
        f.write('camera.XSizeCCD=%d\n' % self.camera['XSizeCCD'])
        f.write('camera.YSizeCCD=%d\n' % self.camera['YSizeCCD'])
        f.write('camera.XBinCCD=%d\n' % self.camera['XBinCCD'])
        f.write('camera.YBinCCD=%d\n' % self.camera['YBinCCD'])
        f.write('sequence=%s\n' % self.sequence['sequence'])
        f.write('sequence.ExposureTime=%d\n' % self.sequence['ExposureTime'])
        f.write('sequence.ExposureRepeat=%d\n' % self.sequence['ExposureRepeat'])
        f.write('sequence.ExposureIncrease=%s\n' % self.sequence['ExposureIncrease'])
        f.write('sequence.FilterSequence=%s\n' % self.sequence['FilterSequence'])
        f.write('sequence.offset=%.1f\n' % self.sequence['offset'])
        f.write('object.ObjectName=%s\n' % self.object['ObjectName'])
        f.write('object.RA=%.5f\n' % self.object['RA'])
        f.write('object.Dec=%.6f\n' % self.object['Dec'])
        f.write('constraints.MoonDistance.Min=%d\n' % self.constraints['MoonDistance.Min'])
        f.write('constraints.SolHeight.Max=%.1f\n' % self.constraints['SolHeight.Max'])
        f.write('constraints.AirmassTarget.Max=%.1f\n' % self.constraints['AirmassTarget.Max'])
        f.write('constraints.AltTarget.Min=%.1f\n' % self.constraints['AltTarget.Min'])
        f.write('file=%s\n' % self.file)
        f.write('duration=%d\n' % self.duration)
        f.flush()
        f.close()

    def fromfile(self):
        """
        loads the data from a save file
        """
        f = open(self.filename, 'rt')
#        lines = f.readlines()
        f.close()
        #TODO: processing of lines
        
    def transfer(self):
        '''
        uploads the files to stella for the submission tool
        
        ssh operator@ciruelo
        rsync -e ssh -r -v  -t sro@stella:/stella/home/www/uploads .

        Dann SCP1 (SCP2):
        SCP1 uploads/weingrill/submit/m_67_bvr_se.xml
        kopiert nach Teneriffa, wifsip
        
        oder FCP1
        findet alle aktuellen (max. 24h alten) dateien und kopiert nach wifsip.
        
        Dann weiter wie im Wiki.
        
        script-update (aus .save .xml files machen ohne GUI):
        
        cd uploads
        ./recreate.sh <liste-von-save-files>
        java -Djava.ext.dirs=/usr/share/java/ -cp /z/operator/java stella.jview.JTargetMaker\$Recreate proposal.create $i
        generiert *hier in upload* die xml files.
        mv *.xml weingrill/submit/,
        dann SCP1, bzw. FCP1
        
        Achtung: Naechstes rsync holt wieder die daten von stella.aip.de
        
        liste-von-save-files: ascii ala:
        
        weingrill/save/ngc_1647_bv_se_20131017.submit
        
        d.h. relativer pfad von uploads aus.

        '''
        from subprocess import call
        import time
        import os
        
        source = self.filename
        target='sro@stella:/stella/home/www/uploads/weingrill/save/'
        time.sleep(1) # otherwise submit.jnlp gets confused
        print 'scp %s %s' % (source, target)
        call(['/usr/bin/scp', source, target])
        print os.path.dirname(source)
        _, filename = os.path.split(source)
        #call(['/usr/bin/ssh', 'operator@ciruelo', '"ls ~/bin/autosubmit.sh %s"' % filename])
        call(['/usr/bin/ssh', 'operator@ciruelo', 'bin/autosubmit.sh %s' % filename])
        
    def tycho(self):
        '''
        plot tycho stars upon fov
        '''
        from tycho import tycho
        
        tycho(self.object['RA'], 
              self.object['Dec'],
              fov = 0.2*1.4142, 
              grid=False, 
              background=False, 
              show=False)


