#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2018 Leibniz-Insitute for Astrophysics Potsdam (AIP)"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__ = "2018-10-09"

from datasource import DataSource
from subprocess import call
import os

def getobjid(objid, targetdir=os.getcwd(), imcopy=False, listonly=False):
    wifsip = DataSource(host='pera', database='stella', user='stella')
    query = "SELECT path, filename "\
            " FROM science "\
            " WHERE objid = '%s'" % objid

    path, filename = wifsip.query(query)[0]
    path = path.replace('stella/wifsip/reduced/','/stella/home/stella/wifsip/reduced/')
    # /stella/home/stella/wifsip/reduced/20140725/
    print("%s/%s" % (path, filename))

    scpsource = os.path.join(path,filename + '.fitz')
    if not listonly:
        call(['scp', 'sro@pina.aip.de:' + scpsource, targetdir])
        source = os.path.join(targetdir, filename+'.fitz[1]')
        target = os.path.join(targetdir, filename+'.fits')
        if imcopy:
            call(['/home/jwe/bin/imcopy', source, target])


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='get WiFSIP frames from pera')
    parser.add_argument('-l', '--list', action='store_true',
                        help='list only, no frames will be downloaded')
    parser.add_argument('-i', '--imcopy', action='store_true',
                        help='convert fitzfile to fits file using imcopy')
    parser.add_argument('objid',
                        help='object ID')
    args = parser.parse_args()

    getobjid(args.objid, imcopy=args.imcopy, listonly=args.list)
