'''
Created on Nov 7, 2013

@author: jwe
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    from opencluster import OpenCluster
    import time

    ngc2422 = OpenCluster(objectname='NGC 2422', uname='NGC 2422 BV', obsmode='bvsl20')           
    print ngc2422.object['RA'], ngc2422.object['Dec']
    ngc2422_subframes = ngc2422.plan_wifsip(nfields=5)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)
    
    for sf in ngc2422_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        time.sleep(1) # otherwise submit.jnlp gets confused
        sf.tofile('/work2/jwe/NGC2422')
        sf.transfer()
    ngc2422.plot(ax)
    ngc2422.tycho()
    time.sleep(1) # otherwise submit.jnlp gets confused
    ngc2422.tofile('/work2/jwe/NGC2422')
    ngc2422.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()
