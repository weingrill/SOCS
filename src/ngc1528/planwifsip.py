'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
Updated May 20, 2016

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_rot(transfer=False):
    import datetime

    ngc1528rot = OpenCluster(objectname='NGC 1528', uname='NGC 1528 rot', obsmode='rot')
    
    ngc1528rot.title = 'SOCS'
    ngc1528rot.abstract = 'Photometric monitoring of open stellar clusters'
    ngc1528rot_subframes = ngc1528rot.plan_wifsip(nfields=4)
    ngc1528rot_group = ClusterGroup(ngc1528rot)
    ngc1528rot_group.startdate =  datetime.datetime(2016, 7, 16) 
    ngc1528rot_group.enddate =  datetime.datetime(2017, 4, 14)
    
    for sf in ngc1528rot_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        ngc1528rot_group.add_daughter(sf.uname)
        if transfer: 
            sf.transfer()
    ngc1528rot_group.tofile(config.projectpath)
    if transfer: 
        ngc1528rot_group.transfer()

def do_cmd(transfer=False):
    ngc1528cmd = OpenCluster(objectname='NGC 1528', uname='NGC 1528', obsmode='BVR')
    ngc1528cmd.title = 'SOCS'
    ngc1528cmd.abstract = 'Photometric monitoring of open stellar clusters'

    ngc1528cmd_subframes = ngc1528cmd.plan_wifsip(nfields=5)
    ngc1528cmd_group = ClusterGroup(ngc1528cmd)
    
    for sf in ngc1528cmd_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        ngc1528cmd_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    ngc1528cmd_group.tofile(config.projectpath)
    if transfer: 
        ngc1528cmd_group.transfer()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='NGC 1528 WiFSIP schedule')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    if args.rot: do_rot(transfer=args.transfer)
    if args.cmd: do_cmd(transfer=args.transfer)
