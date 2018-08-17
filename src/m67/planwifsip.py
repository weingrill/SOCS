#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
from clustergroup import ClusterGroup

def do_uvby(transfer=False, plot=False):
    m67uvby = OpenCluster(objectname='M 67', 
                          uname='M 67 uvby', 
                          obsmode='uvby')
    m67uvby.title = 'SOCS'
    m67uvby.abstract = 'Photometric monitoring of open stellar clusters'

    m67uvby_subframes = m67uvby.plan_wifsip(nfields=5)
    m67uvby_group = ClusterGroup(m67uvby)
    
    for sf in m67uvby_subframes:
        print(sf.uname)
        sf.tofile(config.projectpath)
        m67uvby_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    m67uvby_group.tofile(config.projectpath)
    if transfer: 
        m67uvby_group.transfer()
 
    
def do_hby(transfer=False):
    m67hby = OpenCluster(objectname='M 67', 
                          uname='M 67 Hby', 
                          obsmode='Hby')
    m67hby.title = 'SOCS'
    m67hby.abstract = 'Photometric monitoring of open stellar clusters'

    m67hby_subframes = m67hby.plan_wifsip(nfields=5)
    m67hby_group = ClusterGroup(m67hby)
    
    for sf in m67hby_subframes:
        print(sf.uname)
        sf.tofile(config.projectpath)
        m67hby_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    m67hby_group.tofile(config.projectpath)
    if transfer: 
        m67hby_group.transfer()
 
def do_rot(transfer=False):

    m67rot = OpenCluster(objectname='M 67', uname='M 67 rot', obsmode='rot')
    m67rot.title = 'SOCS'
    m67rot.abstract = 'Photometric monitoring of open stellar clusters'

    m67rot_subframes = m67rot.plan_wifsip(nfields=4)
    m67rot_group = ClusterGroup(m67rot)
    
    for sf in m67rot_subframes:
        print(sf.uname)
        sf.tofile(config.projectpath)
        m67rot_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    m67rot_group.tofile(config.projectpath)
    if transfer: 
        m67rot_group.transfer()
 
def do_cmd(transfer=False):
    m67cmd = OpenCluster(objectname='M 67', uname='M 67 BVR', obsmode='BVR')
    m67cmd.title = 'SOCS'
    m67cmd.abstract = 'Photometric monitoring of open stellar clusters'

    m67cmd_subframes = m67cmd.plan_wifsip(nfields=5)
    m67cmd_group = ClusterGroup(m67cmd)
    
    for sf in m67cmd_subframes:
        print(sf.uname)
        sf.tofile(config.projectpath)
        m67cmd_group.add_daughter(sf.uname)
        if transfer: sf.transfer()
    m67cmd_group.tofile(config.projectpath)
    if transfer: 
        m67cmd_group.transfer()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='M 67 WiFSIP schedule')
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
    
