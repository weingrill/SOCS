#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 11, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import config
from diffphot_nopca import DiffPhotometryNoPCA

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Differential photometry')
    parser.add_argument('-l', '--load',   action='store_true', help='loads the objids from database')
    parser.add_argument('-b', '--build',  action='store_true', help='build the photometry matrix')
    parser.add_argument('-r', '--reduce', action='store_true', help='reduce the photometry matrix')
    parser.add_argument('-c', '--clean',  action='store_true', help='clean the matrix and perform PCA')
    parser.add_argument('-s', '--save',   action='store_true', help='save lightcurves')
    parser.add_argument('-m', '--make',   action='store_true', help='make plots of lightcurves')
    parser.add_argument('--twosigma',   action='store_true', help='two sigma removal of epochs')
    parser.add_argument('field', help='field to process')

    args = parser.parse_args()
    
#    fields = ['M 67 rot NE','M 67 rot NW','M 67 rot SE','M 67 rot SW']
#    for field in fields[:1]:
    diffphot = DiffPhotometryNoPCA(args.field, 
                                   config.datapath, 
                                   config.lightcurvespath, 
                                   config.plotpath)
    if args.load:   diffphot.load_objids()
    if args.build:  diffphot.build_photmatrix()
    if args.reduce: diffphot.reduce()
    if args.clean:  diffphot.clean(twosigma = args.twosigma)
    if args.save:   diffphot.save_lightcurves()
    if args.make:   diffphot.make_lightcurves(show=False)
