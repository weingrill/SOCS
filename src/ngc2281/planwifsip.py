'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_rot(transfer=False):
    from opencluster import OpenCluster

    ngc2281 = OpenCluster(objectname='NGC 2281', 
                          uname='NGC2281 rot', 
                          obsmode='rot')
    ngc2281.title = 'SOCS'
    ngc2281.abstract = 'Photometric monitoring of open stellar clusters'
    ngc2281_subframes = ngc2281.plan_wifsip(nfields=4)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    for sf in ngc2281_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        sf.tofile('/work2/jwe/ngc2281')
        #sf.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()

def do_bvi(transfer=False):
    from opencluster import OpenCluster

    ngc2281 = OpenCluster(objectname='NGC 2281', 
                          uname='NGC 2281 BVI', 
                          obsmode='bvisl')
    ngc2281.title = 'SOCS'
    ngc2281.abstract = 'Photometric monitoring of open stellar clusters'
    ngc2281_subframes = ngc2281.plan_wifsip(nfields=5)
    ngc2281.mode['timeout'] = 5030000.0 # redoing two fields in BVI
    for sf in ngc2281_subframes:
        print sf.uname
        sf.tofile('/work2/jwe/ngc2281')
        sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    #do_rot(transfer=False)
    do_bvi(transfer=True)
