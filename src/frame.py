#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 29, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
from datasource import DataSource


class Frame(object):
    '''
    Frame class to inference the frames on WiFSIP database and filestorage
    '''

    def __init__(self, objid):
        '''
        Constructor
        '''
        table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
        query = """SELECT * FROM frames WHERE objid = '%s';""" % objid
        self._data = table.query(query)[0]

        self.columns = ['star', 'mag_isocor', 'magerr_isocor', 'mag_auto', 'magerr_auto',
                        'flux_auto', 'fluxerr_auto', 'flux_max', 'background',
                        'xwin_image', 'ywin_image',
                        'alphawin_j2000', 'deltawin_j2000',
                        'class_star', 'flags', 'fwhm_image', 'ellipticity']

        params = {'columns': ', '.join(self.columns),
                  'objid': objid}
        query = """SELECT %(columns)s FROM phot WHERE objid = '%(objid)s';""" % params
        self._stars = table.query(query)

        dtypes = ['i4', 'f2', 'f2', 'f2', 'f2',
                  'f4', 'f4', 'f4', 'f4',
                  'f4', 'f4',
                  'f8', 'f8',
                  'f2', 'i4', 'f4', 'f2']
        arraydata = []
        for star in self._stars:
            arraydata.append(tuple(star))
        self.stars = np.array(arraydata, dtype=zip(self.columns, dtypes))

        # self.stars = Stars(objid)
        # print self._data.keys()
        # self.exposure = self._data['expt']

    def __getattr__(self, name):
        '''
        convert the table columns to class properties
        '''
        return self._data[name]

    def download(self, targetdirectory='.'):
        """downloads the file from pina"""
        from subprocess import call
        sourcepath = '/stella/home/stella/wifsip/reduced/' + self.objid[:8] + '/'
        filename = 'science' + self.objid[:14] + 'EXP' + self.objid[-4:] + '.fitz'
        print(filename, sourcepath)
        source = ''.join(['sro@pina.aip.de:', sourcepath, filename, '.fitz '])
        # /stella/home/stella/wifsip/reduced/20130716/science20130716A-0006EXP0011
        call(['scp', source, targetdirectory])
        self.fitzfile = targetdirectory + filename

    def convertfile(self):
        """convert the compressed fitzfile to a fitsfile"""
        from subprocess import call

        self.fitsfile = self.fitzfile[:-1] + 's'
        call(['/home/jwe/bin/imcopy', self.fitzfile, self.fitsfile])

    def __str__(self, *args, **kwargs):
        lines = []
        for column in self._data.keys():
            lines.append('%-14s %s' % (column, self.__getattr__(column)))
        return '\n'.join(lines)


if __name__ == '__main__':
    frame1 = Frame('20130602A-0003-0016')
    print(frame1.expt, frame1.objid)
    print(frame1.stars['alphawin_j2000'])
    print(frame1.stars[0])
    # print frame1
    frame1.download('/work2/jwe/stella')
