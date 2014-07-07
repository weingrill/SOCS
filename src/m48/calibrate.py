'''
Created on Jun 2, 2014

@author: jwe

Tasks:

* load B and V frames from wifsip db
* load tycho-2 ref stars
* match stars
* 

'''

class Calibrate(object):
    '''
    classdocs
    '''


    def __init__(self, filtercol='V'):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        if filtercol in ('B','V'):
            self.filtercol = filtercol
        else:
            raise(ValueError)
    
    def getframes(self, refframe):
        query = "SELECT object FROM frames WHERE objid='%s';" % refframe
        obj = self.wifsip.query(query)[0][0]
        print refframe,'-->',obj
        query = """SELECT objid
         FROM frames
         WHERE object LIKE '%s'
         AND filter like '%s';""" % (obj, self.filtercol)
        
        result = self.wifsip.query(query)
        self.frames = [r[0] for r in result]
        

    def getrefframes(self):
        '''
        choose the frames with the highest number of matched stars for each 
        field and filter 
        '''
        
        query = """SELECT object,max(matched)
            FROM frames
            WHERE object LIKE 'M 48 rot%%'
            AND filter = '%s'
            GROUP BY object
            ORDER BY object;""" % self.filtercol
        result = self.wifsip.query(query)
        
        for r in result:
            print r[0],'\t',r[1]
        
        maxmatched = [{'object':r[0], 'matched':r[1]} for r in result]
        
        self.ref = []
        for m in maxmatched:
            query = """SELECT objid
                FROM frames
                WHERE object LIKE '%s'
                AND filter like '%s'
                AND matched = %d;""" % \
                (m['object'], self.filtercol, m['matched'])
            self.ref.append(self.wifsip.query(query)[0][0])

    def corrframes(self, refframe=''):
        from numpy import array, mean, average, std
        
        query = """CREATE VIEW phot1 AS 
            SELECT * FROM phot WHERE objid='%s'
            ORDER BY phot.mag_auto LIMIT 500;""" % (refframe)
        self.wifsip.execute(query)
        
        for frame in self.frames:
            query = """SELECT phot.mag_auto, phot1.mag_auto
            FROM phot, phot1
            WHERE phot.objid='%s'
            AND circle(phot.coord,0) <@ circle(phot1.coord,0.15/3600.);
            """ % frame
            
            result = self.wifsip.query(query)
            
            omag = array([r[0] for r in result])
            cmag = array([r[1] for r in result])
            wts = 1./(2.512*cmag)
            mstd = std(omag-cmag)
            
            try:
                corr = average(omag-cmag, weights=wts)
            except ZeroDivisionError:
                corr = mean(omag-cmag)
            
            print '%s: %.3f %.3f (%d)' % (frame, corr, mstd, len(omag))
            self.updateframe(frame, corr)
                
        self.wifsip.dropview('phot1')
        
    def updateframe(self, frame, corr):
        from numpy import isnan
        if isnan(corr):
            query = """UPDATE frames SET corr = NULL
                    WHERE objid = '%s';""" % frame
        else:
            query = """UPDATE frames
                        SET corr = %f
                        WHERE objid = '%s';""" % (corr, frame)
        self.wifsip.execute(query)
        
    def resetframes(self):
        query = """UPDATE frames
        SET corr=NULL
        WHERE object LIKE 'M 48 rot%%'
        AND filter like '%s';""" % self.filtercol
        self.wifsip.execute(query)
        
if __name__ == '__main__':
    calv = Calibrate(filtercol='V')
    calv.resetframes()
    calv.getrefframes()
    for r in calv.ref: print r
    
    for ref in calv.ref:
        print 'reference:', ref
        calv.getframes(ref)
        calv.corrframes(ref)
        
#     calb = Calibrate(filtercol='B')
#     calb.resetframes()
#     calb.getrefframes()
#     for r in calb.ref: print r
#     
#     for ref in calb.ref:
#         print 'reference:', ref
#         calb.getframes(ref)
#         calb.corrframes(ref)