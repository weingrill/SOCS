#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 9, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

from opencluster import OpenCluster
from clustergroup import ClusterGroup
class PlanWiFSIP(object):
    '''
    classdocs
    '''


    def __init__(self, objectname, nfields=5, transfer=False):
        '''
        Constructor
        '''
        self.transfer = transfer
        self.objectname = objectname
        self.nfields = nfields
        self.projectpath = '/work2/jwe/'+objectname.replace(' ','')+'/'
        
    def do_uvby(self):
        obs_uvby = OpenCluster(self.objectname, 
                              obsmode='uvby')
        obs_uvby.title = 'SOCS'
        obs_uvby.abstract = 'Photometric monitoring of open stellar clusters'
    
        obs_uvby_subframes = obs_uvby.plan_wifsip(self.nfields)
        obs_uvby_group = ClusterGroup(obs_uvby)
        
        for sf in obs_uvby_subframes:
            print sf.uname
            sf.tofile(self.projectpath)
            obs_uvby_group.add_daughter(sf.uname)
            if self.transfer: sf.transfer()
        obs_uvby_group.tofile(self.projectpath)
        if self.transfer: 
            obs_uvby_group.transfer()
     
        
    def do_hby(self):
        obs_hby = OpenCluster(self.objectname, 
                               obsmode='Hby')
        obs_hby.title = 'SOCS'
        obs_hby.abstract = 'Photometric monitoring of open stellar clusters'
    
        obs_hby_subframes = obs_hby.plan_wifsip(self.nfields)
        obs_hby_group = ClusterGroup(obs_hby)
        
        for sf in obs_hby_subframes:
            print sf.uname
            sf.tofile(self.projectpath)
            obs_hby_group.add_daughter(sf.uname)
            if self.transfer: sf.transfer()
        obs_hby_group.tofile(self.projectpath)
        if self.transfer: 
            obs_hby_group.transfer()
     
    def do_rot(self):
    
        obs_rot = OpenCluster(self.objectname, 
                               obsmode='rot')
        obs_rot.title = 'SOCS'
        obs_rot.abstract = 'Photometric monitoring of open stellar clusters'
        
        # we do not need a central field
        rotfields = self.nfields
        if self.nfields==5:
            rotfields = 4
        
        obs_rot_subframes = obs_rot.plan_wifsip(rotfields)
        obs_rot_group = ClusterGroup(obs_rot)
        
        for sf in obs_rot_subframes:
            print sf.uname
            sf.tofile(self.projectpath)
            obs_rot_group.add_daughter(sf.uname)
            if self.transfer: sf.transfer()
        obs_rot_group.tofile(self.projectpath)
        if self.transfer: 
            obs_rot_group.transfer()
     
    def do_cmd(self):
        obs_cmd = OpenCluster(self.objectname, obsmode='BVR')
        obs_cmd.title = 'SOCS'
        obs_cmd.abstract = 'Photometric monitoring of open stellar clusters'
    
        obs_cmd_subframes = obs_cmd.plan_wifsip(self.nfields)
        obs_cmd_group = ClusterGroup(obs_cmd)
        
        for sf in obs_cmd_subframes:
            print sf.uname
            sf.tofile(self.projectpath)
            obs_cmd_group.add_daughter(sf.uname)
            if self.transfer: sf.transfer()
        obs_cmd_group.tofile(self.projectpath)
        if self.transfer: 
            obs_cmd_group.transfer()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='WiFSIP observation submission')
    parser.add_argument('objectname', 'name of the cluster to observe')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('-uvby', action='store_true', help='uvby observations')
    parser.add_argument('-Hby', action='store_true', help='hahb observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    pw = PlanWiFSIP(args.objectname)
    if args.Hby: pw.do_hby(transfer=args.transfer)
    if args.uvby: pw.do_uvby(transfer=args.transfer)
    if args.rot: pw.do_rot(transfer=args.transfer)
    if args.cmd: pw.do_cmd(transfer=args.transfer)
    
        