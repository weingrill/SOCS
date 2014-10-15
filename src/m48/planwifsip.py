#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

import config
from opencluster import OpenCluster

def do_rot(transfer=False):
    m48 = OpenCluster(objectname='M 48', 
                          uname='M 48 rot', 
                          obsmode='rot')
    m48.title = 'SOCS'
    m48.abstract = 'Photometric monitoring of open stellar clusters'
    m48_subframes = m48.plan_wifsip(nfields=4)

    for sf in m48_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()

def do_cmd(transfer=False):
    m48 = OpenCluster(objectname='M 48', 
                          uname='M 48 BVR', 
                          obsmode='BVR')
    m48.title = 'SOCS'
    m48.abstract = 'Photometric monitoring of open stellar clusters'
    m48_subframes = m48.plan_wifsip(nfields=5)
    for sf in m48_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()
        
if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description='M48 WiFSIP schedule')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    if args.rot: do_rot(transfer=args.transfer)
    if args.cmd: do_cmd(transfer=args.transfer)
