'''
Created on Oct 23, 2013

@author: jwe <jweingrill@aip.de>
'''
from datasource import DataSource

class SuperRef(object):
    def __init__(self, refframe):
        self.ref = {}
        self.refframe = refframe
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        self.offset = {}
        
    def load_master(self):
        """load the master reference frame and populate the ref table"""
        query = """SELECT id, mag_isocor
                   FROM matched, phot
                   WHERE phot.objid like '%s'
                    AND phot.objid = matched.objid
                    AND phot.star = matched.star
                    AND mag_isocor IS NOT NULL;""" % self.refframe
        data = self.wifsip.query(query)
        for d in data:
            try:
                self.ref[d[0]] = d[1]
            except TypeError:
                print d[0]
    
    def reflist(self):
        """return a list of reference frames"""
        query = """select objid 
                   from frames 
                   where objid like '20101005A-0004-00%';"""  
        data = self.wifsip.query(query)
        return data
        
    def correct(self, mag):
        """correct the offset of the mag dict:
        calculate the magnitude weighted average of mag-ref and 
        return the corrected table"""
        from numpy import array, average, mean
        
        corr = {}
        omag = array([mag[key] for key in set(self.ref.keys()) & set(mag.keys())])
        cmag = array([self.ref[key] for key in set(self.ref.keys()) & set(mag.keys())])
        
        weights = 1./omag
        
        avg = average(omag-cmag, weights=weights)
        #avg = mean(omag-cmag)
        print avg
        for key in set(self.ref.keys()) & set(mag.keys()):
            corr[key] = mag[key]-avg
        print len(corr)
        return corr
        
    def balance(self, reflist):
        import numpy as np
        
        # populate avg
        avg = {}
        n = len(reflist)
        for key in self.ref.keys():
            avg[key] = np.empty(n)
            avg[key][:] = np.nan 
        
        for frame in reflist:
            i = reflist.index(frame)
            query = """SELECT id, mag_isocor
                       FROM matched, phot
                       WHERE phot.objid like '%s'
                        AND phot.objid = matched.objid
                        AND phot.star = matched.star
                        AND mag_isocor<=16.5
                        AND mag_isocor>=12.5;""" % frame
            data = self.wifsip.query(query)
            mag = {}
            for d in data:
                mag[d[0]] = d[1]

            # correct the average offset of the list    
            corr = self.correct(mag) 
            #print frame[0]
            #print corr
            
            for key in corr.keys():
                #print key
                try:
                    avg[key][i] = corr[key]
                except ValueError:
                    print key,'has no value'
                except KeyError:
                    print key,' not found'
                    
        self.ref = {}
        self.std = {}
        for key in avg.keys():
            i = np.where(np.isfinite(avg[key]))
            self.ref[key] = np.mean(avg[key][i])
            self.std[key] = np.std(avg[key][i])
            #print key, np.mean(avg[key][i]), np.std(avg[key][i])
    
    def tofile(self, filename):
        """save superref stars to file"""
        from numpy import isfinite
        f = open(filename, 'wt')
        for key in self.ref.keys():
            if isfinite(self.ref[key]):
                f.write('%s\t%f\t%f\n' % (key,self.ref[key],self.std[key]))
        f.close()
    
if __name__ == '__main__':
    sr = SuperRef('20101005A-0004-0006')
    sr.load_master()
    reflist = sr.reflist()
    sr.balance(reflist)
    sr.tofile('/home/jwe/20101005A-0004-0006.tsv') 
