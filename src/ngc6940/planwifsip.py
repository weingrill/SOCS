#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013
Updated May 20, 2016

@author: jwe <jweingrill@aip.de>

start 2016-03-30 jd = 2457477.5
end   2016-12-15 jd = 2457752.5

25.05.2016 18:30:31 [WARN] ObservingHelper…:841 MoonHeight not valid {Max=null}
25.05.2016 18:30:31 [WARN] ObservingHelper…:841 MoonPhase not valid {Max=null}
25.05.2016 18:30:32 [WARN] DOMTarget.parse…:280 No image type specified

'''
import config  # @UnresolvedImport
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_rot(transfer=False):
    import datetime

    ngc6940rot = OpenCluster(objectname='NGC 6940', uname='NGC 6940 rot', obsmode='rot')
    ngc6940rot.priority = 0.5
    ngc6940rot.mode['pernight'] = 6

    ngc6940rot.title = 'SOCS'
    ngc6940rot.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940rot_subframes = ngc6940rot.plan_wifsip(nfields=4)
    ngc6940rot_group = ClusterGroup(ngc6940rot)
    ngc6940rot_group.startdate =  datetime.datetime(2017, 4, 15) 
    ngc6940rot_group.enddate =  datetime.datetime(2017, 12, 29)
    
    for sf in ngc6940rot_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        ngc6940rot_group.add_daughter(sf.uname)
        if transfer: 
            sf.transfer()
    ngc6940rot_group.tofile(config.projectpath)
    if transfer: 
        ngc6940rot_group.transfer()

def do_bvr(transfer=False):
    ngc6940bvr = OpenCluster(objectname='NGC 6940', 
                          obsmode='BVR')
    ngc6940bvr.startdate = '2018-06-06'
    ngc6940bvr.enddate = '2018-07-05'
    ngc6940bvr_subframes = ngc6940bvr.plan_wifsip(nfields=5)
    for sf in ngc6940bvr_subframes:
        print sf.uname
        sf.tofile(config.projectpath+'submit/')
        if transfer:
            sf.transfer()

def do_uvby(transfer=False):
    ngc6940uvby = OpenCluster(objectname='NGC 6940', 
                          obsmode='uvby')
    ngc6940uvby.startdate = '2018-06-06'
    ngc6940uvby.enddate = '2018-07-05'
    ngc6940uvby_subframes = ngc6940uvby.plan_wifsip(nfields=5)
    for sf in ngc6940uvby_subframes:
        print sf.uname
        sf.tofile(config.projectpath+'submit/')
        if transfer:
            sf.transfer()


        
if __name__ == '__main__':
    #do_rot(transfer=True)
    #do_bvr(transfer=True)
    do_uvby(transfer=True)
