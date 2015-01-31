#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 9, 2014

@author: jwe
'''

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
    parser = argparse.ArgumentParser(description='NGC 6633 analysis')
    parser.add_argument('--calibrate', action='store_true', help='create table')
    parser.add_argument('--calibrate2', action='store_true', help='create table')
    parser.add_argument('--create', action='store_true', help='create table')
    parser.add_argument('--clear', action='store_true', help='clear table')
    parser.add_argument('--getframes', action='store_true', help='get frames')
    parser.add_argument('--sigmas', action='store_true', help='plot cmd')
    parser.add_argument('--bv', action='store_true', help='update B-V')
    parser.add_argument('filter', default='V', type=str, help='filter color to process')

    args = parser.parse_args()

    if args.calibrate: calibrate('NGC 6633 BVI', args.filter)
    if args.calibrate2:
        from calibrate2 import Calibrate2
        cal = Calibrate2('NGC 6633 rot NW', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot NE', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot SW', filtername=args.filter)
        cal.grid()
        cal = Calibrate2('NGC 6633 rot SE', filtername=args.filter)
        cal.grid()
        
    if args.create or args.clear or args.getframes or args.sigmas or args.bv:
        from photometry import Photometry
        phot = Photometry(objname='NGC 6633 BVI', filtercol=args.filter, dbname='ngc6633')
    if args.create: phot.createtable()
    if args.clear: phot.cleartable()
    if args.getframes: phot.getframes()
    if args.sigmas: phot.update_sigmas()
    if args.bv: phot.update_bv()