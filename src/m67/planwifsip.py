#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster
def do_uvby(transfer=False, plot=False):
    m67 = OpenCluster(objectname='M 67', 
                          uname='M 67 uvby', 
                          obsmode='uvby')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'
    m67.tofile(config.projectpath)
    if transfer: m67.transfer()

    
def do_hahb(transfer=False):
    m67 = OpenCluster(objectname='M 67', 
                          uname='M 67 hahb', 
                          obsmode='hahb')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'
    m67.tofile(config.projectpath)
    if transfer: m67.transfer()

def do_rot(transfer=False):

    m67 = OpenCluster(objectname='M 67', uname='M 67 rot', obsmode='rot')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'

    m67_subframes = m67.plan_wifsip(nfields=4)
    for sf in m67_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()

def do_cmd(transfer=False):
    m67 = OpenCluster(objectname='M 67', uname='M 67 BVR', obsmode='BVR')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'

    m67_subframes = m67.plan_wifsip(nfields=5)
    for sf in m67_subframes:
        print sf.uname
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='M 67 WiFSIP schedule')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('-uvby', action='store_true', help='uvby observations')
    parser.add_argument('-hahb', action='store_true', help='hahb observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    if args.hahb: do_hahb(transfer=args.transfer)
    if args.uvby: do_uvby(transfer=args.transfer)
    if args.rot: do_rot(transfer=args.transfer)
    if args.cmd: do_cmd(transfer=args.transfer)
    
