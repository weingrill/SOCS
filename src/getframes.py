#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 3, 2013

@author: jwe <jweingrill@aip.de>

get files from M 48 BVI
'''
#TODO: argparse
def getframes(obj, targetdir='/work2/jwe/stella/wifsip/', filtercol='V',
              fwhm=3.0, background=500):
    from datasource import DataSource
    from subprocess import call

    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    query = """SELECT path, filename
             FROM frames, science
             WHERE (filter LIKE '%s')
              AND OBJECT LIKE '%s%%'
              AND frames.objid = science.objid
              AND fwhm_image < %f
              AND backgrnd < %f
              order by frames.objid""" % (filtercol, obj, fwhm, background)
    tab = wifsip.query(query)
    for path, filename in tab:
        print path+'/'+filename
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get frames from pina')
    parser.add_argument('-o', '--object', dest='obj', 
                        help='name of object to fetch' )
    parser.add_argument('-f', '--filter', dest='filtercol', default='V',
                        help='filter color usually V')
    parser.add_argument('--fwhm', type=float, default=3.0,
                        help='filter limit')
    parser.add_argument('--background', type=float, default=500.0,
                        help='background limit')
    parser.add_argument('targetdir', 
                        help='target directory for images')
    args = parser.parse_args()
    
    print args
    getframes(args.obj, targetdir=args.targetdir, filtercol=args.filtercol,
              fwhm=args.fwhm, background=args.background)