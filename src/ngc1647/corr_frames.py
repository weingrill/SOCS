'''
Created on Oct 23, 2013

@author: jwe <jweingrill@aip.de>
'''

class CorrFrames(object):

    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.frames = []
        self.ref = {}
        self.std = {}
        self.corr = {}
        self.wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        
    def load_frames(self):
        query = """select objid
                   from frames
                    where filter like 'rp'
                    and object like 'NGC 1647 field 1'
                    and expt = 62.5
                    and good
                    order by objid"""
        self.frames = self.wifsip.query(query)

    def load_ref(self, filename):
        f = open(filename, 'rt')
        lines = f.readlines()
        f.close()
        for l in lines:
            key, mag, std = l.split('\t')
            self.ref[key] = float(mag)
            self.std[key] = float(std)
    
    def load_stars(self, objid):
        query = """SELECT id, mag_isocor
                   FROM matched, phot
                   WHERE phot.objid like '%s'
                    AND phot.objid = matched.objid
                    AND phot.star = matched.star;""" % objid
        data = self.wifsip.query(query)
        stars = {}
        for d in data:
            stars[d[0]] = d[1]
        return stars
        
    def update_frame(self, objid, corr):
        query = """update frames
                   SET corr = %f
                   WHERE objid like '%s';""" % (corr,objid)
        self.wifsip.execute(query)
        
    def correct(self):
        from numpy import array, average, mean, isfinite
        n = len(self.frames)
        for objid in self.frames:
            i = self.frames.index(objid)
            mag = self.load_stars(objid)
            omag = array([mag[key] for key in set(self.ref.keys()) & set(mag.keys())])
            cmag = array([self.ref[key] for key in set(self.ref.keys()) & set(mag.keys())])
            wts  = 1./array([self.std[key] for key in set(self.std.keys()) & set(mag.keys())])
            try:
                corr = average(omag-cmag, weights=wts)
            except ZeroDivisionError:
                corr = mean(omag-cmag)
            if isfinite(corr):
                print i, n, objid[0], corr
                self.corr['key'] = corr
                self.update_frame(objid[0], corr)

if __name__ == '__main__':
    cf = CorrFrames()
    cf.load_frames()
    cf.load_ref('/home/jwe/20101005A-0004-0006.tsv')
    cf.correct()