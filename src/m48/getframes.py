#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 3, 2013

@author: jwe <jweingrill@aip.de>

get files from M 48 BVI
'''
#TODO: argparse
def getframes(obj, targetdir='/work2/jwe/stella/wifsip/m48'):
    from datasource import DataSource
    from subprocess import call

    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    query = """SELECT path, filename
             FROM frames, science
             WHERE (filter LIKE 'V')
              AND OBJECT LIKE '%s%%'
              AND frames.objid = science.objid
              AND fwhm_image < 3
              AND backgrnd < 500
              order by frames.objid""" % obj
    tab = wifsip.query(query)
    for path, filename in tab:
        print path+'/'+filename
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])

def getrotframes(targetdir='/work2/jwe/stella/wifsip/m48/rot'):
    from datasource import DataSource
    from subprocess import call

    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    query = """SELECT path, filename
             FROM frames, science
             WHERE object LIKE 'M 48 rot%%'
              AND frames.objid = science.objid
              order by frames.objid"""
    tab = wifsip.query(query)
    for path, filename in tab:
        print path+'/'+filename
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='M48 analysis')
    parser.add_argument('--clear', action='store_true', help='clear periods')
    parser.add_argument('--object', dest='obj' )
    getframes()