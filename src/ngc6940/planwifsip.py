'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_rot(transfer=False):
    from opencluster import OpenCluster

    ngc6940 = OpenCluster(objectname='NGC 6940', uname='NGC 6940 rot', obsmode='rot')
    ngc6940.title = 'SOCS'
    ngc6940.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940_subframes = ngc6940.plan_wifsip(nfields=4)

    
    for sf in ngc6940_subframes:
        print sf.uname, sf.duration
        sf.tofile('/work2/jwe/NGC6940')
        if transfer:
            sf.transfer()

def do_bvi(transfer=False):
    from opencluster import OpenCluster

    ngc6940 = OpenCluster(objectname='NGC6940', uname='NGC6940 BVI', obsmode='BVI')
    ngc6940.title = 'SOCS'
    ngc6940.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940_subframes = ngc6940.plan_wifsip(nfields=5)
    for sf in ngc6940_subframes:
        print sf.uname, sf.duration
        sf.tofile('/work2/jwe/NGC6940')
        if transfer:
            sf.transfer()

def do_bvr(transfer=False):
    from opencluster import OpenCluster

    ngc6940 = OpenCluster(objectname='NGC6940', uname='NGC6940 BVR', obsmode='BVR')
    ngc6940.title = 'SOCS'
    ngc6940.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940_subframes = ngc6940.plan_wifsip(nfields=5)
    for sf in ngc6940_subframes:
        print sf.uname, sf.duration
        
        sf.tofile('/work2/jwe/NGC6940')
        if transfer:
            sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    do_rot(transfer=True)
    #do_bvi()
    do_bvr(transfer=True)
