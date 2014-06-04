'''
Created on Mar 19, 2014

@author: jwe
'''

class FrameCal(object):
    '''
    classdocs
    '''


    def __init__(self, filter='V'):
        '''
        Constructor
        '''
        from datasource import DataSource

        self.filter = filter
        self.frames = []
        self.slopes = {}
        self.intercepts = {}
        self.r_values = {}
        self.p_values = {}
        self.std_errs = {} 

        self.wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        #self.frameslist()
        #self.calframes()
        
    def frameslist(self):
        '''
        Get the list of frames in a certain filter
        '''
        query = """SELECT objid, stars
        FROM frames
        WHERE object LIKE 'NGC 2281 BVI %%' 
        AND filter LIKE '%s';""" % self.filter
        result = self.wifsip.query(query)
        self.frames = [r[0] for r in result]
        
        
    def calframe(self, frame):
        if frame in self.slopes:
            return (self.slopes[frame], 
                   self.intercepts[frame],  
                   self.r_values[frame],
                   self.p_values[frame],
                   self.std_errs[frame])
        else:
            raise KeyError('objid %s not found' % frame)
            return (1.0, 0.0, 0.0, 0.0, 0.0)

    def updateframe(self, objid):
        print objid
        query = "UPDATE frames SET good=TRUE WHERE objid = '%s';" % objid
        self.wifsip.execute(query)

    def calframes(self, silent=True):
        import numpy as np
        from scipy import stats
        
        if self.filter == 'V':
            mag = 'vmag'
        if self.filter == 'B':
            mag = 'bmag'
            
        for frame in self.frames:
            query = """SELECT phot.mag_auto, %s
                         FROM frames, phot, ngc2281ref
                         WHERE frames.objid like '%s'
                         AND frames.objid=phot.objid
                         AND frames.filter='%s'
                         AND circle(ngc2281ref.coord,3./3600.) @> circle(phot.coord,.0)
                         AND bv>0.0;""" % (mag, frame, self.filter)
            result = self.wifsip.query(query)
            omags = np.array([s[0] for s in result])
            cmags = np.array([s[1] for s in result])
            slope, intercept, r_value, p_value, std_err = stats.linregress(omags, cmags)    
            self.slopes[frame] = slope
            self.intercepts[frame] = intercept
            self.r_values[frame] = r_value
            self.p_values[frame] = p_value
            self.std_errs[frame] = std_err  
            if std_err < 0.02:
                self.updateframe(frame)
            if not silent:
                print '%s\t%.3f\t%6.3f\t%.3f\t%.4f\t%.3f\t%3d' % \
                (frame, slope, intercept, r_value, p_value, 
                 std_err, len(omags))
        print len(self.frames), ' frames calibrated'
                
        
        
if __name__ == '__main__':
    fcv = FrameCal(filter='V')
    fcv.frameslist()
    fcv.calframes(silent=False)
    
    print fcv.calframe('20140226A-0067-0005')
    fcb = FrameCal(filter='B')
    fcb.frameslist()
    fcb.calframes(silent=False)
    