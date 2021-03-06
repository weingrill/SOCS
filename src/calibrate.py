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
import numpy as np
import psycopg2

try:
    import config
except ImportError:
    print("Please define config file first")
    exit()


class Calibrate(object):
    '''
    classdocs
    '''

    def __init__(self, obj='', filtercol='V'):
        '''
        Constructor:
        
        Initialize database
        set filtercolor
        '''
        from datasource import DataSource

        self.wifsip = DataSource(database='wifsip',
                                 host='pina',
                                 user='sro')
        self.obj = obj
        if filtercol in ('B', 'V', 'R', 'I'):
            self.filtercol = filtercol
        else:
            raise (ValueError)
        logging.basicConfig(filename=config.projectpath + 'calibration.log',
                            format='%(asctime)s %(levelname)s %(message)s',
                            level=logging.INFO)
        logging.info('Object: %s', self.obj)
        logging.info('Setting filtercolor=%s', self.filtercol)

    def getframes(self, refframe):
        logging.info('getframes ...')
        query = "SELECT object FROM frames WHERE objid='%s';" % refframe
        obj = self.wifsip.query(query)[0][0]
        print(refframe, '-->', obj)
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
            WHERE object LIKE '%s'
            AND filter = '%s'
            GROUP BY object
            ORDER BY object;""" % (self.obj, self.filtercol)
        result = self.wifsip.query(query)
        logging.info('%d frames' % len(result))
        for r in result:
            print(r[0], '\t', r[1])

        maxmatched = [{'object': r[0], 'matched': r[1]} for r in result]

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
        logging.info('create view phot1')

        try:
            self.wifsip.execute('DROP VIEW phot1;')
        except psycopg2.ProgrammingError:
            logging.info('drop view phot1 failed')
        query = """CREATE VIEW phot1 AS 
            SELECT * 
            FROM phot 
            WHERE objid='%s'
             AND phot.mag_auto>12 and phot.mag_auto<16
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
            if len(result) == 0:
                logging.warn('no data for frame %s' % frame)
            omag = np.array([r[0] for r in result])
            cmag = np.array([r[1] for r in result])
            mstd = np.std(omag - cmag)

            corr = np.mean(omag - cmag)
            s = '%s %6.3f %.3f %3d' % (frame, corr, mstd, len(omag))

            logging.info(s)
            if len(omag) < 100 or mstd > 0.015:
                corr = np.nan
                # print '#'
            else:
                print(s)
            self.updateframe(frame, corr)

        logging.info('drop view phot1')
        self.wifsip.dropview('phot1')

    def updateframe(self, frame, corr):

        if np.isnan(corr):
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
        WHERE object LIKE '%s'
        AND filter like '%s';""" % (self.obj, self.filtercol)
        self.wifsip.execute(query)
