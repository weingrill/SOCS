#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''
import config  # @UnresolvedImport
from opencluster import OpenCluster

def do_rot(cluster, transfer=False):
    ngc6709 = OpenCluster(objectname = cluster, obsmode = 'rot')
    ngc6709.tofile(config.projectpath)
    if transfer:
        ngc6709.transfer()

def do_frot(cluster, transfer=False):
    ngc6709frot = OpenCluster(objectname = cluster, obsmode = 'frot')
    ngc6709frot.tofile(config.projectpath)
    if transfer:
        ngc6709frot.transfer()
 
 
def do_rottest(cluster, transfer=False):
    ngc6709 = OpenCluster(objectname = cluster, obsmode = 'rottest')
    ngc6709.tofile(config.projectpath)
    if transfer:
        ngc6709.transfer()

            
def do_bvr(cluster, transfer=False):

    ngc6709 = OpenCluster(objectname = cluster, obsmode='BVR')
    ngc6709.tofile(config.projectpath)
    if transfer:
        ngc6709.transfer()
        
if __name__ == '__main__':
    #do_rot('NGC 6709', transfer=True)
    do_frot('NGC 6709', transfer=True)
    #do_rottest('NGC 6709', transfer=True)
    #do_bvr('NGC 6709', transfer=True)
