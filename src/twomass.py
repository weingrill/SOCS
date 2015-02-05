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
            values = self._byname(param)[0] 
        if type(param) is tuple:
            values = self._bycoord(param)[0]
        keys = ['twomass', 'raj2000', 'dej2000', 'jmag', 'e_jmag' ,'hmag' , 
                'e_hmag', 'kmag' , 'e_kmag', 'coord'] 
        for key, value in zip(keys,values):
            if key=='coord': # we want to return a tuple
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
        WHERE point(%f,%f) <@ circle(coord,0.0006) LIMIT 1;""" % coord 
        result = self.corot.query(query)
        return result
        
         
    
tm = TwoMass('18363826+0601086')
print tm['coord']
tm = TwoMass((279.159443,6.01901))
print tm['twomass']
print tm.values()