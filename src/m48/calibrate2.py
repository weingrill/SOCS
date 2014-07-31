'''
Created on Jul 22, 2014

@author: jwe
'''

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
        
    def get_matched(self):
        from numpy import zeros
        starnumbers = zeros(len(self.objids))
        for objid in self.objids:
            print objid,
            query = """SELECT id 
        FROM matched
        WHERE matched.objid='%s'
        ORDER by id;""" % objid
            stars = [r[0] for r in self.wifsip.query(query)]
            print len(stars)
            
            for star in stars:
                self.starids.add(star)
                self.stars[objid] = stars
                i = self.objids.index(objid)
                starnumbers[i] = len(stars) 
        
        print self.epochs, 'number of objids'
        print self.numstars, 'unique stars found'
        
#         from matplotlib import pylab
#         pylab.hist(starnumbers, 20)
#         
#         pylab.savefig('/work2/jwe/m48/stars_'+self.field+'.pdf')

    def create(self):
        from numpy import zeros,nan
        self.a = zeros([self.epochs,self.numstars])
        stararr = [s for s in self.starids]
        for objid in self.objids:
            phot = {}
            query = """
            SELECT id, phot.mag_auto
            FROM frames, phot, matched
            WHERE frames.objid = '%s'
            and frames.objid=phot.objid
            and (phot.objid,phot.star) = (matched.objid,matched.star)
            and phot.flags=0
            """ % objid
            result = self.wifsip.query(query)
            for r in result:
                phot[r[0]] = r[1]
            epoch = self.objids.index(objid)
            for starid in stararr:
                star = stararr.index(starid)
                if starid in phot: 
                    if phot[starid]<30.0:
                        self.a[epoch, star] = phot[starid]
                    else:
                        self.a[epoch, star] = nan
                else:
                    self.a[epoch, star] = nan
                    
        #self.a.dump('/work2/jwe/m48/photmatrix')

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


    def clean(self):
        '''
        remove stars and objids with a large number of nans
        '''
        from numpy import isfinite, delete
        
        delvec = []
        newobjids = []
        for i in range(self.epochs):
            objidvec = self.a[i,:]
            if objidvec[isfinite(objidvec)].size<self.epochs/2:
                delvec.append(i)
            else:
                newobjids.append(self.objids[i])
        delete(self.a, delvec, 0)
        self.objids = newobjids
        print 'deleted ',len(delvec),'objids'
        
        delvec = []
        newstarlist = []
        
        for j in range(self.numstars):
            starvec = self.a[:,j]
            if starvec[isfinite(starvec)].size<self.numstars/2:
                delvec.append(i)
            else:
                newstarlist.append(self.starids[i])
        delete(self.a, delvec, 1)
        self.starids = set(newstarlist)
        print 'deleted ',len(delvec),'stars'
            

    
    def calibrate(self):
        from numpy import load, isfinite
        from scipy.stats import nanmean, nanstd
        
        #self.a = load('/work2/jwe/m48/photmatrix')
        #print self.a.shape
#        self.epochs, self.numstars = self.a.shape
        m = nanmean(self.a, axis=0)
        self.mags = m
        # m now contains the mean magnitude for each star
        for i in range(self.epochs):
            objidvec = self.a[i,:]
            log('/work2/jwe/m48/data/objidvec', '%s %d\n' %
                (self.objids[i], objidvec[isfinite(objidvec)].size))
            self.a[i,:] = self.a[i,:] - m
        
        corr = nanmean(self.a,axis=1)
        self.corr = corr
        
        # corr now contains to corrections for each epoch
        starlist = list(self.starids)
        for j in range(self.numstars):
            starvec = self.a[:,j]
            log('/work2/jwe/m48/data/starvec', '%s %d\n' %
                (starlist[j], starvec[isfinite(starvec)].size))
            self.a[:,j] = self.a[:,j] - corr
        
        std = nanstd(self.a, axis=0)
        for starid, mag,s in zip(self.starids,m,std):
            if isfinite(mag) and isfinite(s):
                log('/work2/jwe/m48/data/stars', '%s %.3f %.4f\n' %(starid, mag, s))
        
    def array_toimage(self):
        from numpy import rint, isnan
        from PIL import Image  # @UnresolvedImport
        from functions import scaleto
        
        a1 = self.a

        a1[isnan(a1)] = 0.0 
        simg = scaleto(a1,[0.0,255.0])
    
        simg = rint(simg)
        simg = simg.astype('uint8')
        im = Image.fromarray(simg)
        
        im.save('/work2/jwe/m48/'+self.field+'.png')                    
        
        #TODO: interate
    def setcorr(self):
        from numpy import nan
        for objid in self.objids:
            i = self.objids.index(objid)
            corr = self.corr[i]
            if abs(corr)>0.03:
                corr = nan
            self.updateframe(objid, corr) 
        
if __name__ == '__main__':
    
    fields = ['M 48 rot NE','M 48 rot NW','M 48 rot SE''M 48 rot SW']
    
    
    for field in fields[0]:
        cal = Calibrate2('M 48 rot NE')
        cal.get_matched()
        cal.create()
        cal.clean()
        cal.calibrate()
        #cal.setcorr()
        cal.array_toimage()
    
    exit()
