'''
Created on Apr 12, 2013
Updated May 20, 2016

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_rot(transfer=False):
    import datetime

    ngc2281rot = OpenCluster(objectname='NGC 2281', uname='NGC 2281 rot', obsmode='rot')
    
    ngc2281rot.title = 'SOCS'
    ngc2281rot.abstract = 'Photometric monitoring of open stellar clusters'
    ngc2281rot_subframes = ngc2281rot.plan_wifsip(nfields=4)
    ngc2281rot_group = ClusterGroup(ngc2281rot)
    ngc2281rot_group.startdate =  datetime.datetime(2016, 8, 20) 
    ngc2281rot_group.enddate =  datetime.datetime(2017, 5, 13)
    
    for sf in ngc2281rot_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        ngc2281rot_group.add_daughter(sf.uname)
        if transfer: 
            sf.transfer()
    ngc2281rot_group.tofile(config.projectpath)
    if transfer: 
        ngc2281rot_group.transfer()

if __name__ == '__main__':
    do_rot(transfer=True)
