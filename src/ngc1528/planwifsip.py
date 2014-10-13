'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
from opencluster import OpenCluster

def do_rot(transfer=False):

    ngc1528 = OpenCluster(objectname='NGC 1528', uname='NGC 1528 rot', obsmode='rot')
    ngc1528.title = 'SOCS'
    ngc1528.abstract = 'Photometric monitoring of open stellar clusters'
    ngc1528_subframes = ngc1528.plan_wifsip(nfields=4)

    
    for sf in ngc1528_subframes:
        print sf.uname, sf.duration
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()

def do_cmd(transfer=False):
    ngc1528 = OpenCluster(objectname='NGC 1528', uname='NGC 1528 BVR', obsmode='BVR')
    ngc1528.title = 'SOCS'
    ngc1528.abstract = 'Photometric monitoring of open stellar clusters'
    ngc1528_subframes = ngc1528.plan_wifsip(nfields=5)
    for sf in ngc1528_subframes:
        print sf.uname, sf.duration
        
        sf.tofile(config.projectpath)
        if transfer: sf.transfer()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='NGC1528 WiFSIP schedule')
    parser.add_argument('-rot', action='store_true', help='rot observations')
    parser.add_argument('-cmd', action='store_true', help='CMD observations')
    parser.add_argument('--transfer', action='store_true', 
                        help='transfert files via sftp to stella')

    args = parser.parse_args()
    
    if args.rot: do_rot(transfer=args.transfer)
    if args.cmd: do_cmd(transfer=args.transfer)
