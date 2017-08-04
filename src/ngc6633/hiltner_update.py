#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 3, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
import config
from ngc6633star import NGC6633Star

def sethiltner():
    with open(config.datapath + 'hiltner_match.csv', 'rt') as f:
        lines = f.readlines()
    for l in lines[1:]:
        starid, hiltner = l.rstrip('\n').split(',')
        print starid, hiltner
        star = NGC6633Star(starid)
        star['hiltner'] = hiltner
        
def sethiltnermember():
    members = []
    nonmembers = []
    with open(config.datapath + 'hiltner.txt', 'rt') as f:
        lines = f.readlines()
    for l in lines[1:]:
        record = l.rstrip('\n').split(' ')
        print record[0], record[-1]
        if record[-1]=='X' or record[-1]=='D':
            members.append(int(record[0]))
        else:
            nonmembers.append(int(record[0]))
    print members
    print nonmembers
    
if __name__ == '__main__':
    sethiltnermember()