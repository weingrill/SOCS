#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 5, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''


class Ppmxl(dict):
    '''
    Class to inference the PPMXL table in the wifsip database on pina
    '''

    def __init__(self, param):
        '''
        Constructor:
        initialize the database and query it according the parameter,
        if a string is given, then we assume, it is the name,
        if it is a tuple, then we take them as coordinates
        '''
        from datasource import DataSource

        self.corot = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        if type(param) is str:
            try:
                values = self._byname(param)[0]
            except IndexError:
                values = [None, None, None, None, None, None, None, None, None, None]

        if type(param) is tuple:
            try:
                values = self._bycoord(param)[0]
            except IndexError:
                values = [None, None, None, None, None, None, None, None, None, None]

        keys = ['id', 'ra', 'dec', 'pmra', 'pmde', 'jmag', 'hmag', 'kmag',
                'b1mag', 'b2mag', 'r1mag', 'r2mag', 'imag', 'smags', 'no',
                'fl', 'coord']
        for key, value in zip(keys, values):
            if key == 'coord':  # we want to return a tuple
                vals = value.strip('()').split(',')
                self[key] = (float(vals[0]), float(vals[1]))
            else:
                self[key] = value

    def _byname(self, name):
        """query the table by name"""
        query = """SELECT * 
        FROM ppmxl 
        WHERE id = '%s';""" % name
        result = self.corot.query(query)
        return result

    def _bycoord(self, coord):
        """query the table by coordinate"""
        query = """SELECT * 
        FROM ppmxl 
        WHERE circle(coord,0.0006) @> circle(point(%f,%f),0) LIMIT 1;""" % coord
        result = self.corot.query(query)
        return result


if __name__ == '__main__':
    tm = Ppmxl('PPMXL5492093740337166978')
    print(tm['coord'])
    tm = Ppmxl((283.961511, 6.064574))
    print(tm['id'])
    print(tm.values())
