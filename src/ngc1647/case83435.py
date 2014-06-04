'''
Created on Dec 3, 2013

@author: jwe

phase dispersion minimization taken from
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyTimingDoc/pyPDMDoc/examples.html#carry-out-a-pdm-analysis

'''
class Case83435(object):
    def __init__(self):
        self.id = 'PPMXL2998058814458783435'
        #self.period = 0.16328
        #self.period = 0.1643
        self.period = 0.328628
    
    def lightcurve_fromdb(self):
        """
        extract a single lightcurve from the database
        and return epoch (hjd), magnitude and error
        """ 
        from datasource import DataSource
        import numpy as np
    
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        query = """SELECT frames.hjd, phot.mag_isocor, phot.magerr_isocor
                FROM frames, matched, phot
                WHERE id LIKE '%s'
                AND filter LIKE 'rp'
                AND frames.objid = matched.objid
                AND (phot.objid,phot.star) = (matched.objid,matched.star)
                
                ORDER BY hjd;""" % (self.id)
    # AND hjd>2455473.5 AND hjd<2455477 AND frames.good
                
    # AND frames.good
    #          AND hjd>2455470 AND hjd<2455510
        data = wifsip.query(query)
        wifsip.close()
        self.hjd = np.array([d[0] for d in data])
        self.mag = np.array([d[1] for d in data])
        self.err = np.array([d[2] for d in data])
        self.hjd -= self.hjd[0]
    
    def plot(self):
        import pylab as plt
        plt.scatter(self.hjd, self.mag)
        plt.show()

    def plot_phased(self):
        import pylab as plt
        import numpy as np
        t = self.hjd % self.period
        i = np.argsort(t)
        t = t[i]
        y = self.mag[i]
        plt.scatter(t, y)
        plt.title(self.id)
        plt.ylim(max(y),min(y))
        plt.grid()
        plt.show()
        
        

    def scanner(self):
        # Import PDM module
        from PyAstronomy.pyTiming import pyPDM
        
        # Get Scanner instance
        scanner = pyPDM.Scanner(minVal=0.1, maxVal=1.0, dVal=0.01, mode="period")
        # Print the periods covered by the scanner
        print "Periods: ",
        for period in scanner:
            print period,    
    
    def analysis(self):
        import numpy
        import matplotlib.pylab as mpl
        from PyAstronomy.pyTiming import pyPDM
        
        # Create artificial data with frequency = 3,
        # period = 1/3
        x = self.hjd
        y = self.mag
        
        # Get a ``scanner'', which defines the frequency interval to be checked.
        # Alternatively, also periods could be used instead of frequency.
        S = pyPDM.Scanner(minVal=0.163, maxVal=0.33, dVal=0.00001, mode="period")
        
        # Carry out PDM analysis. Get frequency array
        # (f, note that it is frequency, because the scanner's
        # mode is ``frequency'') and associated Theta statistic (t).
        # Use 10 phase bins and 3 covers (= phase-shifted set of bins).
        P = pyPDM.PyPDM(x, y)
        p1, t1 = P.pdmEquiBinCover(100, 3, S)
        # For comparison, carry out PDM analysis using 10 bins equidistant
        # bins (no covers).
        p2, t2 = P.pdmEquiBin(100, S)
        
        
        # Show the result
        mpl.figure(facecolor='white')
        mpl.title("Result of PDM analysis")
        mpl.xlabel("Period")
        mpl.ylabel("Theta")
        mpl.plot(p1, t1, 'bp-')
        mpl.plot(p2, t2, 'gp-')
        #mpl.legend(["pdmEquiBinCover", "pdmEquiBin"])
        mpl.show()       

    def my_analysis(self):
        import matplotlib.pylab as mpl
        from pdm import pdm
        
        p1, t1 = pdm(self.hjd, self.mag, 0.162, 0.33, 0.00002, 100)
        mpl.figure(facecolor='white')
        mpl.title("Result of PDM analysis")
        mpl.xlabel("Period")
        mpl.ylabel("Theta")
        mpl.plot(p1, t1)
        mpl.show()       
        
        

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')

    var = Case83435()
    var.lightcurve_fromdb()
    var.plot()
    var.plot_phased()
    #var.analysis()
    var.my_analysis()