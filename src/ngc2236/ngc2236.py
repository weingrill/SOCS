#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 21, 2013

@author: Joerg Weingrill <jweingrill@aip.de>

Data reduction Class for NGC2236
'''


class Ngc2236(object):
    def __init__(self):
        """Constructor"""

    
def calibrate(objname, filtercol):
    from calibrate import Calibrate

    calv = Calibrate(obj=objname+' %%', filtercol=filtercol)
    calv.resetframes()
    calv.getrefframes()
    for r in calv.ref: print r
    
    for ref in calv.ref:
        print 'reference:', ref
        calv.getframes(ref)
        calv.corrframes(ref)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NGC 2236 analysis')
    parser.add_argument('--calibrate', action='store_true', help='calibrate the photometry')
    parser.add_argument('--create', action='store_true', help='create table')
    parser.add_argument('--clear', action='store_true', help='clear table')
    parser.add_argument('--getframes', action='store_true', help='get frames')
    parser.add_argument('--sigmas', action='store_true', help='plot cmd')
    parser.add_argument('--bv', action='store_true', help='update B-V')
    parser.add_argument('filter', default='V', type=str, help='filter color to process')

    args = parser.parse_args()

    if args.calibrate: calibrate('NGC 2236 BVI', args.filter)
        
    if args.create or args.clear or args.getframes or args.sigmas or args.bv:
        from photometry import Photometry
        phot = Photometry(objname='NGC 2236 BVI', filtercol=args.filter, dbname='ngc2236')
    if args.create: phot.createtable()
    if args.clear: phot.cleartable()
    if args.getframes: phot.getframes()
    if args.sigmas: phot.update_sigmas()
    if args.bv: phot.update_bv()
