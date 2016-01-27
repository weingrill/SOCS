#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 18, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
def listframes(obj, fields = ['objid']):
    from datasource import DataSource
    
    wifsip = DataSource(host = 'pera', database = 'stella', user = 'stella')
    params = {'object': obj,
              'fieldsstring': ', '.join(fields)}
    
    query = """SELECT %(fieldsstring)s
             FROM frames
             WHERE object LIKE '%(object)s'
             ORDER by %(fieldsstring)s
             """ % params
    
    tab = wifsip.query(query)
    print '#','\t'.join(fields)
    for t in tab:
        print '\t'.join([str(ti) for ti in t])
    print '#',len(tab),'records'
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='list frames in WiFSIP database')
    parser.add_argument('-f', '--fields', dest='fields', default='filter',
                        help='filter, jd, expt, airmass, fwhm_image, zeropnt, ...')
    parser.add_argument('object', help='name of object')
    args = parser.parse_args()
    fields = args.fields.split(',')
    print args,fields
    listframes(args.object, fields)