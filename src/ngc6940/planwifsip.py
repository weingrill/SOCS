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
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_rot(transfer=False):
    import datetime

    ngc6940rot = OpenCluster(objectname='NGC 6940', uname='NGC 6940 rot', obsmode='rot')
    ngc6940rot.priority = 0.3

    ngc6940rot.title = 'SOCS'
    ngc6940rot.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6940rot_subframes = ngc6940rot.plan_wifsip(nfields=4)
    ngc6940rot_group = ClusterGroup(ngc6940rot)
    ngc6940rot_group.startdate =  datetime.datetime(2016, 7, 15) 
    ngc6940rot_group.enddate =  datetime.datetime(2016, 12, 29)
    
    for sf in ngc6940rot_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        ngc6940rot_group.add_daughter(sf.uname)
        if transfer: 
            sf.transfer()
    ngc6940rot_group.tofile(config.projectpath)
    if transfer: 
        ngc6940rot_group.transfer()

        
if __name__ == '__main__':
    do_rot(transfer=True)
