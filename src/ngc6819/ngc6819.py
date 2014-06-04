'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    from opencluster import OpenCluster
    import time

    ngc6819 = OpenCluster(objectname='NGC 6819', uname='NGC6819 BV', obsmode='bvsl20')
    ngc6819.title = 'SOCS'
    ngc6819.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6819_subframes = ngc6819.plan_wifsip(nfields=5)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)
    
    #ngc6819.mode['timeout'] = 5030000.0 # redoing two fields in BVI
    for sf in ngc6819_subframes:
        print sf.uname
        time.sleep(1) # otherwise submit.jnlp gets confused
        print sf.timeout
        sf.tofile('/work2/jwe/NGC6819')
        sf.transfer()
        #sf.plot(ax)
        #sf.tycho()

    print ngc6819.timeout, ngc6819.duration
    #xmin,xmax = plt.xlim()
    #ax.set_xlim(xmax,xmin)
    #plt.show()
