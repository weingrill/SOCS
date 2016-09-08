'''
Created on Oct 3, 2013

@author: jwe <jweingrill@aip.de>
'''

def airmass(altitude):
    """converts altitude to airmass"""
    from numpy import cos, radians

    return 1./cos(radians(90.-altitude))
    
def altitude(airmass):
    """converts airmass to altitude"""
    from numpy import arccos, degrees
    
    return 90.-degrees(arccos(1./airmass))
    

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
        self.constraints = {'MoonDistance.Min': 30,
                            'SolHeight.Max':   -16.0,
                            'AirmassTarget.Max': 2.0,
                            'AltTarget.Min':    30.0,
                            'MoonHeight.Max':   90.0,
                            'MoonPhase.Max':     1.0}
        self.file = ''
        self.filename = ''
        self.fields = 1
        self.obsmode = obsmode
        self._setobsmode(obsmode)
        
        if self.object['RA'] is None or self.object['Dec'] is None:
            self.get_coordiantes()
    
    def _setobsmode(self, obsmode):
        """
        sets sequence and mode parameters according to obsmode
        """
        
        obsmodes = {'BVI': {'ExposureTime':     24.0,
                            'ExposureIncrease': '1,5,25,1,5,25,1,3,9',
                            'FilterSequence':   'B,B,B,V,V,V,I,I,I'
                            },
                    'BVR': {'ExposureTime':     24.0,
                            'ExposureIncrease': '1,5,25,1,5,25,1,5,25',
                            'FilterSequence':   'B,B,B,V,V,V,R,R,R'
                            },
                    
                    'UBVRI': {'ExposureTime':    24.0,
                            'ExposureIncrease': '1,5,25,1,5,25,1,5,25,1,5,25,1,5,25',
                            'FilterSequence':   'U,U,U,B,B,B,V,V,V,R,R,R,I,I,I',
                            'pernight':         1
                            },
                    
                    'uvby': {'ExposureTime':    24.0,
                            'ExposureIncrease': '1,1,1,1,5,5,5,5,25,25,25,25',
                            'FilterSequence':   'u,v,b,y,u,v,b,y,u,v,b,y',
                            'MoonDistance.Min': 15,
                            'pernight':         1
                            },
                    
                    'Hby': {'ExposureTime':     24.0,
                            'ExposureIncrease': '5,25,5,25,5,25,5,25,5,25,5,25',
                            'FilterSequence':   'b,b,y,y,hbw,hbw,hbn,hbn,haw,haw,han,han',
                            'MoonDistance.Min': 15,
                            'pernight':         1
                            },
                    
                    'rot': {'ExposureTime':     24.0,
                            'ExposureIncrease': '1,5,25',
                            'FilterSequence':   'V,V,R',
                            'pernight':         6 
                            }
                            
                    }

        if not (obsmode in obsmodes):
            raise TypeError('%s not a valid observation mode' % obsmode)
        
        # set defaults:
        # zerofraction is the length of one exposure in days
        self.mode['zerofraction'] = 1.0/24.0                                    
        self.mode['impact']       = 1.0
        self.mode['mode']         = 'Clusters'
        self.mode['pernight']     = 2 
        
        # get the observation parameters depending on the observation mode
        obsparams = obsmodes[obsmode]
        
        # set the mandatory keys
        for key in ['ExposureTime', 'ExposureIncrease', 'FilterSequence']:
            self.sequence[key]     = obsparams[key]
        
        # set the mandatory pernight value
        if 'pernight' in obsparams:
            self.mode['pernight'] = obsparams['pernight']
        
        filtersequencelen = len(self.sequence['FilterSequence'].split(','))
        exposureincreaselen = len(self.sequence['ExposureIncrease'].split(','))
        # filter sequence and exposures must be equal
        assert(filtersequencelen == exposureincreaselen)
        # set ExposureRepeat
        self.sequence['ExposureRepeat'] = filtersequencelen
        # period_day is related to pernight
        self.mode['period_day'] = 0.5/self.mode['pernight'] # was 0.25
        # set the Minimum Moon distance
        if 'MoonDistance.Min' in obsparams:
            self.constraints['MoonDistance.Min'] = obsparams['MoonDistance.Min']
            
        if 'AirmassTarget.Max' in obsparams:
            self.constraints['AirmassTarget.Max'] = obsparams['AirmassTarget.Max']
            self.constraints['AltTarget.Min'] = altitude(self.constraints['AirmassTarget.Max'])
        
        if 'AltTarget.Min' in obsparams:
            self.constraints['AltTarget.Min'] = obsparams['AltTarget.Min']
            self.constraints['AirmassTarget.Max'] = airmass(self.constraints['AltTarget.Min'])
        
    
    def plot_ephem(self, obsdate=None):
        import matplotlib.pyplot as plt
        
        import ephem
        import datetime
        from numpy import pi,empty
        
        #from astronomy import airmass
        
        stella = ephem.Observer()
        #stella.lon, stella.lat = '13.104659', '52.404963' # Potsdam
        #stella.lat, stella.lon = 31.9583, -111.59867 # KPNO
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
        from cluster import Cluster
        
        c = Cluster(self.object['ObjectName'])
        
        self.ra_str,  self.dec_str = c.coordinatestring
        self.object['RA'] = c['ra']
        self.object['Dec'] = c['dec']
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
        def sbool(s):
            if s.lower() == 'true': return True
            elif s.lower() == 'false': return False
 
        with open(self.filename, 'rt') as f:
            lines = f.readlines()
        #f.close()
        floatparams = ['priority',
                       'mode.period_day',
                       'mode.zerofraction'
                       'mode.impact',
                       'object.RA',
                       'object.Dec',
                       'constraints.MoonDistance.Min',
                       'constraints.SolHeight.Max',
                       'constraints.AirmassTarget.Max',
                       'constraints.AltTarget.Min']
        intparams = ['mode.pickdelay', 
                     'mode.pernight', 
                     'camera.XOffCCD',
                     'camera.YOffCCD',
                     'camera.XSizeCCD',
                     'camera.YSizeCCD',
                     'camera.XBinCCD',
                     'camera.YBinCCD',
                     'sequence.ExposureTime',
                     'sequence.ExposureRepeat',
                     'duration']
        boolparams = ['withfocus','withacquire','withguiding']
        parameters = {}
        for l in lines:
            key, value = l.split('=')
            if key in floatparams:
                value = float(value)
            elif key in intparams:
                value = int(value)
            elif key in boolparams:
                value = sbool(value)
            parameters[key] = value
           
        self.startdate =                parameters['startdate']
        self.enddate =                  parameters['enddate']
        self.priority =                 parameters['priority']
        self.telescope =                parameters['TELESCOPE']
        self.withfocus =                parameters['withfocus']
        self.withacquire =              parameters['withacquire']
        self.withguiding =              parameters['withguiding']
        self.title =                    parameters['title']
        self.uname =                    parameters['uname']
        self.propid =                   parameters['propid']
        self.abstract =                 parameters['abstract']
        self.pi =                       parameters['pi']
        self.affil =                    parameters['affil']
        self.team =                     parameters['team']
        self.mode =                     parameters['mode']
        self.mode.pickdelay =           parameters['mode.pickdelay']
        self.mode.pernight =            parameters['mode.pernight']
        self.mode.period_day =          parameters['mode.period_day']
        self.mode.zerofraction =        parameters['mode.zerofraction']
        self.mode.impact =              parameters['mode.impact']
        self.camera =                   parameters['camera']
        self.camera.XOffCCD =           parameters['camera.XOffCCD']
        self.camera.YOffCCD =           parameters['camera.YOffCCD']
        self.camera.XSizeCCD =          parameters['camera.XSizeCCD']
        self.camera.YSizeCCD =          parameters['camera.YSizeCCD']
        self.camera.XBinCCD =           parameters['camera.XBinCCD']
        self.camera.YBinCCD =           parameters['camera.YBinCCD']
        self.sequence =                 parameters['sequence']
        self.sequence.ExposureTime =    parameters['sequence.ExposureTime']
        self.sequence.ExposureRepeat =  parameters['sequence.ExposureRepeat']
        self.sequence.ExposureIncrease = parameters['sequence.ExposureIncrease']
        self.sequence.FilterSequence =  parameters['sequence.FilterSequence']
        self.sequence.offset =          parameters['sequence.offset']
        self.object.ObjectName =        parameters['object.ObjectName']
        self.object.RA =                parameters['object.RA']
        self.object.Dec =               parameters['object.Dec']
        self.constraints.MoonDistance.Min=parameters['constraints.MoonDistance.Min']
        self.constraints.SolHeight.Max = parameters['constraints.SolHeight.Max']
        self.constraints.AirmassTarget.Max=parameters['constraints.AirmassTarget.Max']
        self.constraints.AltTarget.Min = parameters['constraints.AltTarget.Min']
        self.file =                     parameters['file']
        self.duration =                 parameters['duration']
        
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
        
        Dann weiter wie im Wiki:
        http://stella.aip.de/groups/activity/wiki/index.php/Adding_targets
        
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

        stella@wifsip:~/stella/master1/testing/targets
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
        print'executing operator@ciruelo:autosubmit.sh %s' % filename
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
