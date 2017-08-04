#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 19, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''


if __name__ == '__main__':
    from gyroage import gyroperiod
    from numpy import linspace
    #bv = linspace(0.0, 1.6, num=50)
    bv = linspace(0.35, 1.6, num=100)
    
    age = 600
    
    P11 = gyroperiod(bv, age)
    P34 = gyroperiod(bv, age, P0=3.4)
    P01 = gyroperiod(bv, age, P0=0.1)
    
    
    
    for k in zip(bv,P01,P11,P34):
        print '%.3f\t%.3f\t%.3f\t%.3f' % k
