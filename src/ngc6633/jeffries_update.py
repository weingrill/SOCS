#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 3, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
import config
from ngc6633star import NGC6633Star

def setjeffries():
    with open(config.datapath + 'jeffries_match.csv', 'rt') as f:
        lines = f.readlines()
    for l in lines[1:]:
        starid, jeffries = l.rstrip('\n').split(',')
        print starid, jeffries.lstrip('J')
        star = NGC6633Star(starid)
        star['jeffries'] = jeffries.lstrip('J')
        
    
if __name__ == '__main__':
    setjeffries()