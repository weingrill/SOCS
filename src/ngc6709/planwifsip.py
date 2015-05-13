#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster

def do_rot(cluster, transfer=False):
    ngc6709 = OpenCluster(objectname = cluster, obsmode = 'rot')
    ngc6709.title = 'SOCS'
    ngc6709.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6709.tofile(config.projectpath)
    if transfer:
        ngc6709.transfer()
            
def do_bvr(cluster, transfer=False):

    ngc6709 = OpenCluster(objectname = cluster, obsmode='BVR')
    ngc6709.title = 'SOCS'
    ngc6709.abstract = 'Photometric monitoring of open stellar clusters'
    ngc6709.tofile(config.projectpath)
    if transfer:
        ngc6709.transfer()
        
if __name__ == '__main__':
    do_rot('NGC 6709', transfer=True)
    do_bvr('NGC 6709', transfer=True)
