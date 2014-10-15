'''
Created on Jul 22, 2014

@author: jwe
'''
import pickle
try:
    import config
except ImportError:
    print "Please define config file first"
    exit()

def log(filename, message):
    f = open(filename, 'a')
    f.write(message)
    f.close()

class Calibrate2(object):
    '''
    classdocs
    '''


    def __init__(self, field, filtername='V'):
        '''
        Constructor
        fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE','M 48 rot SW']
        '''
        from datasource import DataSource
        self.field = field
        self.filtername = filtername
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        self._getobjids()
        self.stars = {}
        self.starids = set()
        self.datafilebase = config.datapath+self.field+'_'+self.filtername+'_'
        self.plotfilebase = config.datapath+self.field+'_'+self.filtername
        print self.field
    
    def _getobjids(self):
        query = """SELECT objid,hjd 
        FROM frames
        WHERE object like '%s%%'
        AND filter='%s'
        ORDER by objid;""" % (self.field, self.filtername)
        result = self.wifsip.query(query)
        self.objids = [o[0] for o in result]
        self.hjds = [o[1] for o in result]
    
    @property
    def epochs(self):
        return len(self.objids)   
    
    @property
    def numstars(self):
        return len(self.starids)
        
    def get_matched(self, verbose=False, store=True):
        from numpy import zeros
        
        print 'getting matched stars ...'
        
        try:
            picklefile = open(self.datafilebase+'starids.pickle', 'rb')
            self.starids = pickle.load(picklefile)
            picklefile.close()
            picklefile = open(self.datafilebase+'stars.pickle', 'rb')
            self.stars = pickle.load(picklefile)
            picklefile.close()
            
        except:
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
            if store:
                picklefile = open(self.datafilebase+'starids.pickle', 'wb')
                pickle.dump(self.starids, picklefile)
                picklefile.close()
                picklefile = open(self.datafilebase+'star.pickle', 'wb')
                pickle.dump(self.stars, picklefile)
                picklefile.close()    
        print self.epochs, 'objids'
        print self.numstars, 'unique stars found'
        
    def create(self, store=False):
        def dot():
            import sys
            sys.stdout.write('.')
        from numpy import zeros,nan,load, save
        
        print 'creating array ...'
        
        self.a = zeros([self.epochs,self.numstars])
        try:
            self.a = load(self.datafilebase+'photmatrix.npy')
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
                and phot.flags<8;
                """ % objid
                result = self.wifsip.query(query)
                for r in result:
                    phot[r[0]] = r[1]
                epoch = self.objids.index(objid)
                for starid in stararr:
                    if stararr.index(starid) % 100 == 99: 
                        dot()
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
            if store: save(self.datafilebase+'photmatrix.npy', self.a)
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

    def clip(self, sigma=3.0, store=False, verbose=True):
        """
        perform sigma clipping on each lightcurve
        """
        from numpy import save,load, nan
        from scipy.stats import nanstd, nanmean
        try:
            self.a = load(self.datafilebase+'clippedmatrix.npy')
        except IOError:
            print 'Clipping ...',self.numstars,'stars'
            
            #if a star  has datapoints exceeding sigma times std: make them nan 
            for j in range(self.numstars):
                std = nanstd(self.a[:,j])
                m = nanmean(self.a[:,j])
                if verbose:
                    print self.a[:,j], sigma*std
                    #print len(self.a[:,j] > sigma*nanstd(self.a[:,j])), len(self.a[:,j])
                self.a[abs(self.a[:,j]-m) > sigma*std,j] = nan
                
            if store: save(self.datafilebase+'clippedmatrix.npy', self.a)    

    def clean(self, verbose=False, store=False):
        '''
        remove stars and objids with a large number of nans
        '''
        from numpy import isfinite, delete, save,load
        try:
            self.a = load(self.datafilebase+'cleanedmatrix.npy')
        except IOError:
            print 'Cleaning ...',self.epochs,'epochs'

            delvec0 = []
            newobjids = []
            
            #if an epoch shows less than 10% of the stars: remove it
            for i in range(self.epochs):
                objidvec = self.a[i,:]
                if verbose: print objidvec[isfinite(objidvec)].size,self.numstars
                if objidvec[isfinite(objidvec)].size < self.numstars/10:
                    delvec0.append(i)
                else:
                    newobjids.append(self.objids[i])
            #print delvec0
            
            
            print 'Cleaning ...',self.numstars,'stars'
            
            delvec1 = []
            newstarlist = []
            
            #if a star appears in less than 10% of the epochs: remove it 
            for j in range(self.numstars):
                starvec = self.a[:,j]
                if verbose: 
                    print starvec[isfinite(starvec)].size, self.epochs
                if starvec[isfinite(starvec)].size < self.epochs/10:
                    delvec1.append(j)
                else:
                    newstarlist.append(list(self.starids)[j])
            #print delvec1
            self.a = delete(self.a, delvec0, 0)
            self.a = delete(self.a, delvec1, 1)
            self.objids = newobjids
            print 'deleted ',len(delvec0),'objids'
            self.starids = set(newstarlist)
            print 'deleted ',len(delvec1),'stars'
            if store: save(self.datafilebase+'cleanedmatrix.npy', self.a)    

    
    def calibrate(self):
        from numpy import isfinite
        from scipy.stats import nanmean, nanstd
        
        print 'Calibration ...'
        m = nanmean(self.a, axis=0)
        self.mags = m
        # m now contains the mean magnitude for each star
        for i in range(self.epochs):
            objidvec = self.a[i,:]
            log(self.datafilebase+'objidvec', '%s %d\n' %
                (self.objids[i], objidvec[isfinite(objidvec)].size))
            self.a[i,:] = self.a[i,:] - m
        
        corr = nanmean(self.a,axis=1)
        self.corr = corr
        
        # corr now contains to corrections for each epoch
        starlist = list(self.starids)
        for j in range(self.numstars):
            starvec = self.a[:,j]
            log(self.datafilebase+'starvec', '%s %d\n' %
                (starlist[j], starvec[isfinite(starvec)].size))
            self.a[:,j] = self.a[:,j] - corr
        
        std = nanstd(self.a, axis=0)
        for starid, mag,s in zip(self.starids,m,std):
            if isfinite(mag) and isfinite(s):
                log(self.datafilebase+'stars', '%s %.3f %.4f\n' % 
                    (starid, mag, s))
        
        self.storeplot(m, std)
        
    def storeplot(self, magnitudes, errors):
        from matplotlib import pyplot
        
        ax = pyplot.subplot(111) 
        ax.set_yscale('log') 
        ax.scatter(magnitudes, errors)
        pyplot.xlabel('V (mag)')
        pyplot.ylabel('err (mag)')
        pyplot.xlim(min(magnitudes), max(magnitudes))
        pyplot.ylim(min(errors), max(errors))
        pyplot.title(self.field)
        pyplot.savefig(self.plotfilebase+'.pdf')
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

    def grid(self):
        self.get_matched()
        self.create(store=True)
        self.clean()
        self.clip(sigma=1.0)
        self.calibrate()
        self.setcorr()
        self.array_toimage()
        
