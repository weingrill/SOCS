'''
Created on Apr 12, 2013
Updated May 20, 2016

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup
import datetime

def do_rot(transfer=False):

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

def do_phot(transfer=False):

    ngc2281phot = OpenCluster(objectname = 'NGC 2281', obsmode='UBVRI')
    ngc2281phot.title = 'SOCS'
    ngc2281phot.startdate = datetime.datetime(2017, 4, 26) 
    ngc2281phot.enddate = datetime.datetime(2017, 5, 13) 
    
    ngc2281phot.abstract = 'Photometric monitoring of open stellar clusters'
    ngc2281phot_subframes = ngc2281phot.plan_wifsip(nfields=5)
    for sf in ngc2281phot_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        if transfer: 
            sf.transfer()
    


if __name__ == '__main__':
    #do_rot(transfer=True)
    do_phot(transfer=False)
