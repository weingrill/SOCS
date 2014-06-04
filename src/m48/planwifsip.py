'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_rot(transfer=False):
    from opencluster import OpenCluster

    m48 = OpenCluster(objectname='M 48', 
                          uname='M 48 rot', 
                          obsmode='rot')
    m48.title = 'SOCS'
    m48.abstract = 'Photometric monitoring of open stellar clusters'
    m48_subframes = m48.plan_wifsip(nfields=4)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    for sf in m48_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        sf.tofile('/work2/jwe/m48')
        if transfer:
            sf.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()

def do_bvi(transfer=False):
    from opencluster import OpenCluster

    m48 = OpenCluster(objectname='M 48', 
                          uname='M 48 BVI', 
                          obsmode='bvisl')
    m48.title = 'SOCS'
    m48.abstract = 'Photometric monitoring of open stellar clusters'
    m48_subframes = m48.plan_wifsip(nfields=5)
    for sf in m48_subframes:
        print sf.uname
        sf.tofile('/work2/jwe/m48')
        if transfer:
            sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    do_rot(transfer=True)
    #do_bvi(transfer=False)
