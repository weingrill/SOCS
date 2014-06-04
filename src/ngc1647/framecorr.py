#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on Feb 15, 2013

@author: jweingrill@aip.de

corrects frames with the weighted average of all photometric measurements
'''
class FrameCorr(object):
    def __init__(self):
        """Constructor: make database connection"""
        import psycopg2
        self.wifsip = psycopg2.connect(database='wifsip', user='sro', host='pina.aip.de') 
        self.cur = self.wifsip.cursor()  

    def getframes(self, objects):
        """get frames"""
        query = "SELECT objid, filter FROM frames WHERE object LIKE '"+objects+"%' ORDER BY objid;"
        self.cur.execute(query)
        data = self.cur.fetchall()
        return data

    def corrframe(self, objid, flt):
        """correct every frame for variations"""
        from numpy import mean,average,array, isnan, nan
        
        #transform filters
        f = {'U':'u_', 'B':'b_', 'V':'v_', 'R':'r_', 'I':'i_'}
        if flt in f: 
            flt=f[flt]
        # select stars in frame
        query = "select stars.id, stars."+flt+", stars."+flt+"err, phot.mag_isocor, phot.magerr_isocor"+\
                " from stars, matched, phot"+\
                " where matched.objid like '"+objid+"'"+\
                " and matched.id = stars.id"+\
                " and (matched.objid,matched.star)=(phot.objid,phot.star)"+\
                " and stars."+flt+"err > 0.0;"
        self.cur.execute(query)
        data = self.cur.fetchall()
        #id = [datum[0] for datum in data]
        cmag = array([datum[1] for datum in data])
        cmagerr = array([datum[2] for datum in data])
        omag = array([datum[3] for datum in data])
        #omagerr = array([datum[4] for datum in data])
        # calculate weighted average of the stars in the frame
        # correction offset = weighted average - median magnitude
        corr = nan
        # determine their median magnitude (get from starsview)
        try:
            corr = average(omag-cmag,weights=1./cmagerr)
        except ZeroDivisionError:
            corr = mean(omag-cmag)
        except TypeError:
            corr = nan
        finally:
            # store the corrected magnitude in frame
            print '%s (%s): %6.3f' % (objid,flt, corr)
            if not isnan(corr):
                query ="update frames set corr = %f where objid = '%s';" % (corr, objid)
                self.cur.execute(query)
                self.wifsip.commit()
    
    def close(self):
        self.wifsip.close()
    
if __name__ == '__main__':
    # for each frame in frames do
    import sys
    if len(sys.argv)>1:
        objects = sys.argv[1]
    else:
        objects = '%'
    fc = FrameCorr()
    data = fc.getframes(objects)
    for datum in data:
        #print '%s: %s' % (datum[0], datum[1])
        sys.stdout.write('[%d/%d] ' % (data.index(datum)+1,len(data)))
        fc.corrframe(datum[0], datum[1])
    fc.close()