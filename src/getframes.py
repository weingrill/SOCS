#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on May 3, 2013

@author: Jörg Weingrill <jweingrill@aip.de>

get files from M 48 BVI
'''
def getframes(obj, targetdir='/work2/jwe/stella/wifsip/', filtercol='V',
              conditions=None, imcopy=False, listonly=False):
    from datasource import DataSource
    from subprocess import call
    
    params = {'object': obj,
              'filtercol': filtercol}
    
    wifsip = DataSource(host = 'pera', database = 'stella', user = 'stella')
    query = """SELECT path, filename, fwhm, backgrnd, airmass
             FROM frames, science
             WHERE (filter LIKE '%(filtercol)s')
              AND object LIKE '%(object)s'
              AND frames.objid = science.objid""" % params
    if 'fwhm' in conditions:
        params['fwhm'] = conditions['fwhm']
        query += '\nAND fwhm_image < %(fwhm)f ' % params
    if 'background' in conditions:
        params['background'] = conditions['background']
        query += '\nAND backgrnd < %(background)f ' % params
    if 'airmass' in conditions:
        params['airmass'] = conditions['airmass']
        query += '\nAND airmass < %(airmass)f ' % params
    if 'expt' in conditions:
        params['expt'] = conditions['expt']
        query += '\nAND expt = %(expt)f ' % params  
                
    query = query + '\nORDER by frames.objid;'
    tab = wifsip.query(query)
    print len(tab),'files'
    for path, filename, fwhm, backgrnd, airmass in tab:
        print "%s/%s %.1f %.1f %.2f" % (path, filename, fwhm, backgrnd, airmass)
        # if listonly is set, we break the loop here
        if listonly: continue
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])
        source = '%s/%s.fitz[1]' % (targetdir,filename)
        target = '%s/%s.fits' % (targetdir,filename)
        if imcopy:
            call(['/home/jwe/bin/imcopy', source, target])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get WiFSIP frames from pera')
    parser.add_argument('-o', '--object', dest='obj', 
                        help='name of object to fetch' )
    parser.add_argument('-f', '--filter', dest='filtercol', default='V',
                        help='filter color usually V')
    parser.add_argument('-w', '--fwhm', type=float, help='fwhm limit')
    parser.add_argument('-b', '--background', type=float, help='background limit')
    parser.add_argument('-a', '--airmass', type=float, help='airmass limit')
    parser.add_argument('-e', '--expt', type=float, help='exposuretime')
    
    parser.add_argument('-l', '--list', action='store_true',
                        help='list only, no frames will be downloaded')
    parser.add_argument('--imcopy', action='store_true',
                        help='convert fitzfile to fits file using imcopy')
    parser.add_argument('targetdir', 
                        help='target directory for images')
    args = parser.parse_args()
    
    main_conditions = {}
    if args.fwhm:       main_conditions['fwhm']       = args.fwhm
    if args.background: main_conditions['background'] = args.background
    if args.airmass:    main_conditions['airmass']    = args.airmass
    if args.expt:       main_conditions['expt']       = args.expt
    #print args
    getframes(args.obj, targetdir=args.targetdir, filtercol=args.filtercol,
              conditions = main_conditions, imcopy=args.imcopy, listonly=args.list)
