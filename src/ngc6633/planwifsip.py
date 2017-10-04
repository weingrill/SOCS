'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

import config

def do_rot(transfer=False):
    from opencluster import OpenCluster

    ngc6633 = OpenCluster(objectname='NGC 6633', 
                          uname='NGC 6633 rot', 
                          obsmode='rot')
    ngc6633.title = 'SOCS'
    ngc6633.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6633_subframes = ngc6633.plan_wifsip(nfields=4)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    for sf in ngc6633_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        sf.tofile(config.projectpath+'submit/')
        if transfer:
            sf.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    #plt.show()
    plt.savefig(config.plotpath+'ngc6633submit.pdf')
    plt.close()

def do_bvi(transfer=False):
    from opencluster import OpenCluster

    ngc6633 = OpenCluster(objectname='NGC 6633', 
                          obsmode='BVI')
    ngc6633.startdate = '2017-08-25'
    ngc6633.enddate = '2017-09-08'
    ngc6633_subframes = ngc6633.plan_wifsip(nfields=5)
    #ngc6633.mode['timeout'] = 5030000.0 # redoing two fields in BVI
    for sf in ngc6633_subframes:
        print sf.uname
        sf.tofile(config.projectpath+'submit/')
        if transfer:
            sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    #do_rot(transfer=True)
    do_bvi(transfer=False)
