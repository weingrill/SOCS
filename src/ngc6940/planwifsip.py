'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_rot():
    from opencluster import OpenCluster
    import time

    ngc6940 = OpenCluster(objectname='NGC 6940', uname='NGC 6940 rot', obsmode='rot')
    ngc6940.title = 'SOCS'
    ngc6940.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940_subframes = ngc6940.plan_wifsip(nfields=4)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    for sf in ngc6940_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        time.sleep(1) # otherwise submit.jnlp gets confused
        sf.tofile('/work2/jwe/NGC6940')
        #sf.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()

def do_bvi():
    from opencluster import OpenCluster
    import time

    ngc6940 = OpenCluster(objectname='NGC6940', uname='NGC6940 BVI', obsmode='BVI')
    ngc6940.title = 'SOCS'
    ngc6940.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940_subframes = ngc6940.plan_wifsip(nfields=5)
    ngc6940.mode['timeout'] = 5030000.0 # redoing two fields in BVI
    for sf in ngc6940_subframes:
        print sf.uname
        time.sleep(1) # otherwise submit.jnlp gets confused
        sf.tofile('/work2/jwe/NGC6940')
        sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    #do_rot()
    do_bvi()
