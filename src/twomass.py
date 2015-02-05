#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 5, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class TwoMass(dict):
    '''
    Class to inference the TwoMass table in the corot database on pina
    '''

    def __init__(self, param):
        '''
        Constructor:
        initialize the database and query it according the parameter,
        if a string is given, then we assume, it is the name,
        if it is a tuple, then we take them as coordinates
        '''
        from datasource import DataSource
        
        self.corot = DataSource(database='corot', user='sro', host='pina.aip.de')
        if type(param) is str:
            try:
                values = self._byname(param)[0] 
            except IndexError:
                values = [None,None,None,None,None,None,None,None,None,None]
                
        if type(param) is tuple:
            try:
                values = self._bycoord(param)[0]
            except IndexError:
                values = [None,None,None,None,None,None,None,None,None,None]
                
        keys = ['twomass', 'raj2000', 'dej2000', 'jmag', 'e_jmag' ,'hmag' , 
                'e_hmag', 'kmag' , 'e_kmag', 'coord'] 
        for key, value in zip(keys,values):
            if key=='coord' and not value is None: # we want to return a tuple
                vals = value.strip('()').split(',')
                self[key] = (float(vals[0]),float(vals[1])) 
            else:
                self[key] = value    
        
    def _byname(self, name):
        """query the table by name"""
        query = """SELECT * 
        FROM twomass 
        WHERE twomass = '%s';""" % name
        result = self.corot.query(query)
        return result
    
    def _bycoord(self, coord):
        """query the table by coordinate"""
        query = """SELECT * 
        FROM twomass 
        WHERE circle(coord,0.0006) @> point(%f,%f) LIMIT 1;""" % coord 
        result = self.corot.query(query)
        return result
        
if __name__ == '__main__':
    tm = TwoMass('18363826+0601086')
    print tm['coord']
    tm = TwoMass((210.159443,6.01901)) #) (97.487622, 6.888813)
    print tm['twomass']
    print tm.values()