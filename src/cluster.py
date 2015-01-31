#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Cluster(dict):
    '''
    Cluster object to access the clusters table on wifsip database
    '''
    def __init__(self, name):
        '''
        Constructor:
        set the name of the cluster and open database
        '''
        from datasource import DataSource
        
        self.name = name
        
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        
    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'clusters';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        query = """SELECT * from clusters where name = '%s'""" % self.name
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        return values

    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM clusters 
        WHERE name='%s';""" % (key, self.name))
        return result[0][0]    

if __name__ == '__main__':
    c = Cluster('NGC 2236')
    print c['ra'],c['dec']
    print c['d'],c['diam']        