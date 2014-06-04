'''
Created on Dec 5, 2013

@author: jwe
'''

def fileexists(filename):
    try:
        f = open(filename, 'rt')
    except IOError:
        return False
    f.close()
    return True 

def loadfromfile(filename):
    f = open(filename, 'rt')
    lines = f.readlines()
    f.close()
    result = [l.rstrip('\n') for l in lines]
    return result

def savetofile(filename, table):
    f = open(filename, 'wt')
    for item in table:
        f.write('%s\n'% str(item))
    f.close()
    
class StarCal(object):
    '''
    classdocs
    '''


    def __init__(self, filtercol='None'):
        '''
        Constructor
        '''
        self.stars = []
        self.objids = []
        self.filter = filtercol
        
        self.get_frames()
        self.get_stars()
    
    def get_frames(self):
        if fileexists('/work2/jwe/NGC2281/'+self.filter+'frames.txt'):
            self.objids = loadfromfile('/work2/jwe/NGC2281/'+self.filter+'frames.txt')
            print len(self.objids), 'frames'
            return
        
        from datasource import DataSource
     
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        query = """SELECT frames.objid 
                    FROM frames
                    WHERE frames.object like 'NGC 2281 BVI %%'
                    AND frames.filter like '%s';""" % self.filter
        result = wifsip.query(query)
        self.objids = [s[0] for s in result]
        savetofile('/work2/jwe/NGC2281/'+self.filter+'frames.txt', self.objids)
        print len(self.objids), 'frames'
        wifsip.close()
    
    def get_stars(self):
        
        if fileexists('/work2/jwe/NGC2281/'+self.filter+'stars.txt'):
            self.stars = loadfromfile('/work2/jwe/NGC2281/'+self.filter+'stars.txt')
            print len(self.stars), 'stars'
            return
        
        from datasource import DataSource
     
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        query = """SELECT distinct matched.id
            FROM frames, phot, matched
            WHERE frames.object like 'NGC 2281 BVI %%'
            AND frames.filter like '%s'
            AND frames.objid = phot.objid
            AND (matched.objid,matched.star) = (phot.objid,phot.star);"""  % self.filter
        result = wifsip.query(query)
        self.stars = [s[0] for s in result]
        savetofile('/work2/jwe/NGC2281/'+self.filter+'stars.txt', self.stars)
        wifsip.close()

        print len(self.stars), 'stars'

    def update_db(self):
        from datasource import DataSource
     
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        query = """INSERT INTO ngc2281stars
            SELECT distinct matched.id
            FROM frames, phot, matched
            WHERE frames.object like 'NGC 2281 BVI %%'
            AND frames.filter like '%s'
            AND frames.objid = phot.objid
            AND (matched.objid,matched.star) = (phot.objid,phot.star)
            AND matched.id NOT IN (SELECT id FROM ngc2281stars);"""  % self.filter
        wifsip.execute(query)
        wifsip.close()

    def buildtable(self):
        """
        builds the table of stars
        """
        import numpy as np
        
        epochs = len(self.objids)
        stars = len(self.stars)
        if fileexists('/work2/jwe/NGC2281/'+self.filter+'array.npy'):
            m = np.load('/work2/jwe/NGC2281/'+self.filter+'array.npy')
        else:
            from datasource import DataSource
            from framecal import FrameCal
            
            fc = FrameCal(self.filter)
            
            m = np.zeros([epochs, stars])
            # objid is specific to a filter so we only need to query the objid
            wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
            for objid in self.objids:
                k = self.objids.index(objid)
                print k, epochs, objid,
                query = """SELECT matched.id, phot.mag_auto, phot.mag_errauto 
                        FROM phot, matched
                        WHERE phot.objid like '%s'
                        AND (matched.objid,matched.star) = (phot.objid,phot.star)
                        AND phot.flags = 0;""" % objid
                result = wifsip.query(query)
                starids = [s[0] for s in result]
                mags = [s[1] for s in result]
                err = [s[2] for s in result]
                
                slope, intercept, _, _, _ = fc.calframe(objid) 
                print len(mags)
                for starid in starids:
                    i = self.stars.index(starid)
                    m[k, i] = mags[starids.index(starid)]*slope + intercept
            np.save('/work2/jwe/NGC2281/'+self.filter+'array.npy', m)
            wifsip.close()
        
        i = np.where(m == 0.0)
        m[i] = np.nan
        from scipy import stats
        # calculate the observed average for the stars
        avg = stats.nanmean(m, axis=0)
        for k in range(epochs):
            print k,epochs,self.objids[k]
            
            # calculate the mean of offsets 
            off = stats.nanmedian(m[k,:]-avg)
            # correct epoch for mean of offsets
            m[k,:] += off
        
        # calculate new corrected means
        avg = stats.nanmean(m, axis=0)
        std = stats.nanstd(m, axis=0)
        for i in range(len(self.stars)):
            print self.stars[i],avg[i],std[i]


    def update_star(self, datasource, id, mag, err, n):
        if self.filter=='V':
            query = """
            UPDATE ngc2281stars
            SET vmag = %f, vmag_err = %f, nv = %d
            WHERE id like '%s';""" % (mag, err, n, id)
        elif self.filter=='B':
            query = """
            UPDATE ngc2281stars
            SET bmag = %f, bmag_err = %f, nb = %d
            WHERE id like '%s';""" % (mag, err, n, id)
        query = query.replace('nan', 'NULL')    
        datasource.execute(query)

    def phot(self):
        """
        builds the table of stars
        """
        import numpy as np
        from scipy import stats
        
        epochs = len(self.objids)
        stars = len(self.stars)
        from datasource import DataSource
        from framecal import FrameCal
        
        fc = FrameCal(self.filter)
        fc.frameslist()
        fc.calframes()
        
        m = np.zeros([epochs, stars])
        # objid is specific to a filter so we only need to query the objid
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        for star in self.stars:
            print star,
            query = """SELECT frames.objid, mag_auto, magerr_auto 
                    FROM frames, phot, matched
                    WHERE matched.id like '%s'
                    AND frames.filter like '%s'
                    AND frames.objid = phot.objid
                    AND (matched.objid,matched.star) = (phot.objid,phot.star)
                    AND phot.flags = 0
                    AND frames.good
                    AND magerr_auto > 0.0;""" % (star,self.filter)
            result = wifsip.query(query)
            objids = ([s[0] for s in result])
            mags = np.array([s[1] for s in result])
            err = np.array([s[2] for s in result])

            #print mags
            #print err
            if len(mags)>1:
                for objid in objids:
                    i = objids.index(objid)
                    try:
                        slope, intercept, _, _, _ = fc.calframe(objid)
                        mags[i] = mags[i]*slope + intercept
                    except KeyError:
                        mags[i] = np.NaN 
                m = stats.nanmean(mags)
                s = stats.nanstd(mags)
                merr = stats.nanmean(err)
                stderr = stats.nanstd(err)

                print '%4d %.3f %.3f %.3f %.3f' % (len(mags),m,s,merr,stderr),
                mags = mags[err<=merr+stderr]
                err = err[err<=merr+stderr]
                avg = np.average(mags, weights=1./err)
                std =  np.sqrt(np.average(abs(mags-avg)**2, weights=1./err))
                #std = np.std(mags)
                
                 
                
                print '%4d %.3f %.3f' % (len(mags), avg, std)
                self.update_star(wifsip, star, avg, std, len(mags))
            else: 
                print 'none' 
            
        wifsip.close()
    
    def calibrate(self):
        from datasource import DataSource
        
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        wifsip.execute('UPDATE ngc2281stars SET bv = Bmag-Vmag;')
        wifsip.close()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    
    print 'V-Filter'
    vstars = StarCal(filtercol='V')
    vstars.update_db()
    vstars.phot()
    print 'B-Filter'
    bstars = StarCal(filtercol='B')
    bstars.update_db()
    bstars.phot()
    
    bstars.calibrate()