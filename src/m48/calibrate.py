#!/usr/bin/python
'''
Created on Jun 2, 2014

@author: jwe

Tasks:

* load B and V frames from wifsip db
* load tycho-2 ref stars
* match stars
* 

'''
import logging

class Calibrate(object):
    '''
    classdocs
    '''


    def __init__(self, filtercol='V'):
        '''
        Constructor:
        
        Initialize database
        set filtercolor
        '''
        from datasource import DataSource
        
        self.wifsip = DataSource(database='wifsip', 
                                 host='pina', 
                                 user='sro')
        if filtercol in ('B','V'):
            self.filtercol = filtercol
        else:
            raise(ValueError)
        logging.basicConfig(filename='/work2/jwe/m48/m48_calibration.log', 
                            format='%(asctime)s %(levelname)s %(message)s',
                            level=logging.INFO)
        logging.info('Setting filtercolor=%s',self.filtercol)
    
    def getframes(self, refframe):
        logging.info('getframes ...')
        query = "SELECT object FROM frames WHERE objid='%s';" % refframe
        obj = self.wifsip.query(query)[0][0]
        print refframe,'-->',obj
        logging.info('%s --> %s' % (refframe, obj))
        query = """SELECT objid
         FROM frames
         WHERE object LIKE '%s'
         AND filter like '%s';""" % (obj, self.filtercol)
        
        result = self.wifsip.query(query)
        self.frames = [r[0] for r in result]
        logging.info('%d frames' % len(self.frames))

    def getrefframes(self):
        '''
        choose the frames with the highest number of matched stars for each 
        field and filter 
        '''
        logging.info('getrefframes ...')
        query = """SELECT object,max(matched)
            FROM frames
            WHERE object LIKE 'M 48 rot%%'
            AND filter = '%s'
            GROUP BY object
            ORDER BY object;""" % self.filtercol
        result = self.wifsip.query(query)
        logging.info('%d frames' % len(result))
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
        logging.info('%d frames' % len(self.ref))
        
    def corrframes(self, refframe=''):
        from numpy import array, mean, average, std, nan
        logging.info('create view phot1')
        query = """CREATE VIEW phot1 AS 
            SELECT * 
            FROM phot 
            WHERE objid='%s'
             AND phot.mag_auto>10 and phot.mag_auto<16
             AND flags=0
            ORDER BY phot.mag_auto;""" % (refframe)
        self.wifsip.execute(query)
        
        for frame in self.frames:
            query = """SELECT phot.mag_auto, phot1.mag_auto
            FROM phot, phot1
            WHERE phot.objid='%s'
            AND circle(phot.coord,0) <@ circle(phot1.coord,0.15/3600.);
            """ % frame
            
            result = self.wifsip.query(query)
            if len(result)==0:
                logging.warn('no data for frame %s' % frame)
            omag = array([r[0] for r in result])
            cmag = array([r[1] for r in result])
            wts = 1./(2.512*cmag)
            mstd = std(omag-cmag)
            
            try:
                #corr = average(omag-cmag, weights=wts)
                corr = mean(omag-cmag)
            except ZeroDivisionError:
                corr = mean(omag-cmag)
            if len(omag)<250 or mstd>0.5:
                corr = nan
            s = '%s: %.3f %.3f (%d)' % (frame, corr, mstd, len(omag))
            print s
            logging.info(s)
            self.updateframe(frame, corr)
                
        logging.info('drop view phot1')
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
        logging.info('reset frames')
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