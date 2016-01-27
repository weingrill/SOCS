#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_uvby(transfer=False, plot=False):
    ngc2422uvby = OpenCluster(objectname='NGC 2422', 
                          uname='NGC 2422 uvby', 
                          obsmode='uvby')
    ngc2422uvby.title = 'SOCS'
    ngc2422uvby.abstract = 'Photometric monitoring of open stellar clusters'

    ngc2422uvby_subframes = ngc2422uvby.plan_wifsip(nfields=5)
    ngc2422uvby_group = ClusterGroup(ngc2422uvby)
    
    for sf in ngc2422uvby_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        ngc2422uvby_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    ngc2422uvby_group.tofile(config.projectpath)
    if transfer: 
        ngc2422uvby_group.transfer()
 
    
def do_hby(transfer=False):
    ngc2422hby = OpenCluster(objectname='NGC 2422', 
                          uname='NGC 2422 Hby', 
                          obsmode='Hby')
    ngc2422hby.title = 'SOCS'
    ngc2422hby.abstract = 'Photometric monitoring of open stellar clusters'

    ngc2422hby_subframes = ngc2422hby.plan_wifsip(nfields=5)
    ngc2422hby_group = ClusterGroup(ngc2422hby)
    
    for sf in ngc2422hby_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        ngc2422hby_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    ngc2422hby_group.tofile(config.projectpath)
    if transfer: 
        ngc2422hby_group.transfer()
 
def do_rot(transfer=False):

    ngc2422rot = OpenCluster(objectname='NGC 2422', uname='NGC 2422 rot', obsmode='rot')
    '''
    override priority, so it doesn't overrule ngc2323
    '''
    ngc2422rot.priority = 0.8 
    ngc2422rot.title = 'SOCS'
    ngc2422rot.abstract = 'Photometric monitoring of open stellar clusters'

    ngc2422rot_subframes = ngc2422rot.plan_wifsip(nfields=4)
    ngc2422rot_group = ClusterGroup(ngc2422rot)
    
    for sf in ngc2422rot_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        ngc2422rot_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    ngc2422rot_group.tofile(config.projectpath)
    if transfer: 
        ngc2422rot_group.transfer()
 
def do_cmd(transfer=False):
    ngc2422cmd = OpenCluster(objectname='NGC 2422', uname='NGC 2422 BVR', obsmode='BVR')
    ngc2422cmd.title = 'SOCS'
    ngc2422cmd.abstract = 'Photometric monitoring of open stellar clusters'

    ngc2422cmd_subframes = ngc2422cmd.plan_wifsip(nfields=5)
    ngc2422cmd_group = ClusterGroup(ngc2422cmd)
    
    for sf in ngc2422cmd_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        ngc2422cmd_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    ngc2422cmd_group.tofile(config.projectpath)
    if transfer: 
        ngc2422cmd_group.transfer(config.projectpath)
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='NGC2422 WiFSIP schedule')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('-uvby', action='store_true', help='uvby observations')
    parser.add_argument('-Hby', action='store_true', help='hahb observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    if args.Hby: do_hby(transfer=args.transfer)
    if args.uvby: do_uvby(transfer=args.transfer)
    if args.rot: do_rot(transfer=args.transfer)
    if args.cmd: do_cmd(transfer=args.transfer)
    
