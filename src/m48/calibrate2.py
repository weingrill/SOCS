'''
Created on Jul 22, 2014

@author: jwe
'''

import config

def log(filename, message):
    f = open(filename, 'a')
    f.write(message)
    f.close()

class Calibrate2(object):
    '''
    classdocs
    '''


    def __init__(self, field):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.field = field
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        self._getobjids()
        self.stars = {}
        self.starids = set()
        print self.field
    
    def _getobjids(self):
        query = """SELECT objid,hjd 
        FROM frames
        WHERE object like '%s%%'
        AND filter='V'
        ORDER by objid;""" % self.field
        result = self.wifsip.query(query)
        self.objids = [o[0] for o in result]
        self.hjds = [o[1] for o in result]
    
    @property
    def epochs(self):
        return len(self.objids)   
    
    @property
    def numstars(self):
        return len(self.starids)
        
    def get_matched(self, verbose=False):
        from numpy import zeros,load
        
        print 'getting matched stars ...'
        
        try:
            self.stars = load(config.datapath+self.field+'_stars.npy')
            self.starids = load(config.datapath+self.field+'_starids.npy')
        except IOError:
            starnumbers = zeros(len(self.objids))
            for objid in self.objids:
                if verbose: print objid,
                query = """SELECT id 
                FROM matched
                WHERE matched.objid='%s'
                ORDER by id;""" % objid
                stars = [r[0] for r in self.wifsip.query(query)]
                if verbose: print len(stars)
                
                for star in stars:
                    self.starids.add(star)
                    self.stars[objid] = stars
                    i = self.objids.index(objid)
                    starnumbers[i] = len(stars) 
            
            #save(config.datapath+self.field+'_stars.npy', self.stars)
            #save(config.datapath+self.field+'_starids.npy', self.starids)
        print self.epochs, 'objids'
        print self.numstars, 'unique stars found'
        
    def create(self):
        from numpy import zeros,nan,load, save
        
        print 'creating array ...'
        
        self.a = zeros([self.epochs,self.numstars])
        try:
            self.a = load(config.datapath+self.field+'_photmatrix.npy')
        except IOError:
            stararr = [s for s in self.starids]
            assert(len(stararr)>0)
            
            for objid in self.objids:
                print '%.2f %s' % (100.0*self.objids.index(objid)/len(self.objids), objid),
                phot = {}
                query = """
                SELECT id, phot.mag_auto
                FROM frames, phot, matched
                WHERE frames.objid = '%s'
                and frames.objid=phot.objid
                and (phot.objid,phot.star) = (matched.objid,matched.star)
                and phot.flags=0;
                """ % objid
                result = self.wifsip.query(query)
                for r in result:
                    phot[r[0]] = r[1]
                epoch = self.objids.index(objid)
                for starid in stararr:
                    if stararr.index(starid) % 100 == 99: 
                        print '.',
                    star = stararr.index(starid)
                    if starid in phot: 
                        if phot[starid]<30.0:
                            self.a[epoch, star] = phot[starid]
                        else:
                            self.a[epoch, star] = nan
                    else:
                        self.a[epoch, star] = nan
                print '.'
                
            print 'saving photometric matrix ...'            
            save(config.datapath+self.field+'_photmatrix.npy', self.a)
        print 'shape of photometric matrix: ', self.a.shape

    def updateframe(self, frame, corr):
        from numpy import isnan
        print frame,corr
        if isnan(corr):
            query = """UPDATE frames SET corr = NULL
                    WHERE objid = '%s';""" % frame
        else:
            query = """UPDATE frames
                        SET corr = %f
                        WHERE objid = '%s';""" % (corr, frame)
        self.wifsip.execute(query)


    def clean(self, verbose=False):
        '''
        remove stars and objids with a large number of nans
        '''
        from numpy import isfinite, delete, save,load
        from scipy.stats import nanstd
        try:
            self.a = load(config.datapath+self.field+'_cleanedmatrix.npy')
        except IOError:
            print 'Cleaning ...',self.epochs,'epochs'

            delvec0 = []
            newobjids = []
            
            #if an epoch shows less than half of the stars: remove it
            for i in range(self.epochs):
                objidvec = self.a[i,:]
                if verbose: print objidvec[isfinite(objidvec)].size,self.numstars
                if objidvec[isfinite(objidvec)].size < self.numstars/2:
                    delvec0.append(i)
                else:
                    newobjids.append(self.objids[i])
            #print delvec0
            self.objids = newobjids
            print 'deleted ',len(delvec0),'objids'
            
            
            print 'Cleaning ...',self.numstars,'stars'
            
            delvec1 = []
            newstarlist = []
            stderr = nanstd(self.a)
            #if a star appears in less than half of the epochs: remove it 
            for j in range(self.numstars):
                starvec = self.a[:,j]
                if verbose: print starvec[isfinite(starvec)].size, self.epochs
                if starvec[isfinite(starvec)].size < self.epochs/2 \
                or nanstd(starvec) > 3.0*stderr:
                    delvec1.append(j)
                else:
                    newstarlist.append(list(self.starids)[j])
            #print delvec1
            self.a = delete(self.a, delvec0, 0)
            self.a = delete(self.a, delvec1, 1)
            self.starids = set(newstarlist)
            print 'deleted ',len(delvec1),'stars'
            save(config.datapath+self.field+'_cleanedmatrix.npy', self.a)    

    
    def calibrate(self):
        from numpy import isfinite
        from scipy.stats import nanmean, nanstd
        
        print 'Calibration ...'
        m = nanmean(self.a, axis=0)
        self.mags = m
        # m now contains the mean magnitude for each star
        for i in range(self.epochs):
            objidvec = self.a[i,:]
            log(config.datapath+self.field+'_objidvec', '%s %d\n' %
                (self.objids[i], objidvec[isfinite(objidvec)].size))
            self.a[i,:] = self.a[i,:] - m
        
        corr = nanmean(self.a,axis=1)
        self.corr = corr
        
        # corr now contains to corrections for each epoch
        starlist = list(self.starids)
        for j in range(self.numstars):
            starvec = self.a[:,j]
            log(config.datapath+self.field+'_starvec', '%s %d\n' %
                (starlist[j], starvec[isfinite(starvec)].size))
            self.a[:,j] = self.a[:,j] - corr
        
        std = nanstd(self.a, axis=0)
        for starid, mag,s in zip(self.starids,m,std):
            if isfinite(mag) and isfinite(s):
                log(config.datapath+self.field+'_stars', '%s %.3f %.4f\n' % 
                    (starid, mag, s))
        
        self.storeplot(m, std)
        #TODO: interate
        
    def storeplot(self, magnitudes, errors):
        from matplotlib import pyplot
        
        pyplot.scatter(magnitudes, errors)
        pyplot.xlabel('V (mag)')
        pyplot.ylabel('err (mag)')
        pyplot.xlim(9.0, 21.0)
        pyplot.ylim(0.001, 0.3)
        pyplot.title(self.field)
        pyplot.savefig(config.plotpath+self.field+'.pdf')
        pyplot.close()
        
        
    def array_toimage(self):
        from numpy import rint, isnan
        from PIL import Image  # @UnresolvedImport
        from functions import scaleto
        
        import pyfits
        
        hdu = pyfits.PrimaryHDU(self.a)
        hdu.writeto(config.plotpath+self.field+'.fits', clobber=True)
        a1 = self.a
        
        a1[isnan(a1)] = 0.0
        simg = scaleto(a1,[0.0,255.0])
    
        simg = rint(simg)
        simg = simg.astype('uint8')
        im = Image.fromarray(simg)
        
        im.save(config.plotpath+self.field+'.png')                    
        
        
    def setcorr(self):
        from numpy import nan
        for objid in self.objids:
            i = self.objids.index(objid)
            corr = self.corr[i]
            if abs(corr)>0.03:
                corr = nan
            self.updateframe(objid, corr) 
        
if __name__ == '__main__':
    
    fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE','M 48 rot SW']
    
    
    for field in fields:
        cal = Calibrate2(field)
        cal.get_matched()
        cal.create()
        cal.clean()
        cal.calibrate()
        cal.setcorr()
        cal.array_toimage()
    
    exit()
