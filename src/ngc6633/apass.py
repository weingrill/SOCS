#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Dec 1, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
class Apass(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def fromfile(self, filename=None):
        import numpy as np
        if filename is None:
            filename = config.datapath+'apass.csv'
        self.table = np.genfromtxt(filename, delimiter=',',names=True, dtype=(float, float, float, int, float, float, float, float, int, int, int, 
                                                                              float, float, float, float, float, float, float, float, float, float, float, float))
        
    
    def todatabase(self):
        from datasource import DataSource
        import StringIO
        
        values =''
        
        for t in self.table:
            values = values + 'APASS%(recno)d\t%(RAJ2000).6f\t%(DEJ2000).6f\t%(Bmag).3f\t%(Vmag).3f\t%(rmag).3f\t%(imag).3f\n' % t
            
        wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        f = StringIO.StringIO(values)
        columns = ['starid','ra', 'dec', '"B"', '"V"','"R"','"I"']
        wifsip.cursor.copy_from(f,'referencestars',columns=columns, null='nan')
        wifsip.commit()
        wifsip.close()
        
        
        
apass = Apass()
apass.fromfile()
apass.todatabase()
