#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 29, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np

class Frame(object):
    '''
    Frame class to inference the frames on WiFSIP database and filestorage
    '''
    def __init__(self, objid):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        table = DataSource(database='stella', user='stella', host='pera.aip.de', dictcursor=True)
        query = """SELECT * FROM frames WHERE objid = '%s';""" % objid
        self._data = table.query(query)[0]
        
        self.columns = ['star', 'mag_isocor', 'magerr_isocor', 'mag_auto', 'magerr_auto', 
                   'flux_auto', 'fluxerr_auto', 'flux_max', 'background', 
                   'xwin_image','ywin_image', 
                   'alphawin_j2000',  'deltawin_j2000', 
                   'class_star', 'flags',   'fwhm_image', 'ellipticity']

        params = {'columns':', '.join(self.columns),
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
        self.stars = np.array(arraydata, dtype = zip(self.columns, dtypes))
        
        
        #self.stars = Stars(objid)
        #print self._data.keys()
        #self.exposure = self._data['expt']
        
    def __getattr__(self, name):
        '''
        convert the table columns to class properties
        '''
        return self._data[name]
    
    def download(self, targetdirectory):
        pass
    
    def __str__(self, *args, **kwargs):
        lines = []
        for column, data in zip(self.columns, self._data):
            lines.append('%s\t%s' % (column, data))
        return '\n'.join(lines)
    
        
if __name__ == '__main__':
    frame1 = Frame('20130602A-0003-0016')
    print frame1.expt, frame1.objid
    print frame1.stars['alphawin_j2000']
    print frame1.stars[0]
    print frame1
    