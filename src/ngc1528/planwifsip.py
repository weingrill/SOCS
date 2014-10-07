'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_rot(transfer=False):
    from opencluster import OpenCluster

    ngc1528 = OpenCluster(objectname='NGC 1528', uname='NGC 1528 rot', obsmode='rot')
    ngc1528.title = 'SOCS'
    ngc1528.abstract = 'Photometric monitoring of open stellar clusters'
    ngc1528_subframes = ngc1528.plan_wifsip(nfields=4)

    
    for sf in ngc1528_subframes:
        print sf.uname, sf.duration
        sf.tofile('/work1/jwe/SOCS/NGC1528')
        if transfer:
            sf.transfer()

def do_bvr(transfer=False):
    from opencluster import OpenCluster

    ngc1528 = OpenCluster(objectname='NGC 1528', uname='NGC 1528 BVR', obsmode='BVR')
    ngc1528.title = 'SOCS'
    ngc1528.abstract = 'Photometric monitoring of open stellar clusters'
    ngc1528_subframes = ngc1528.plan_wifsip(nfields=5)
    for sf in ngc1528_subframes:
        print sf.uname, sf.duration
        
        sf.tofile('/work1/jwe/SOCS/NGC1528')
        if transfer:
            sf.transfer()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    do_rot(transfer=True)
    #do_bvi()
    do_bvr(transfer=True)
