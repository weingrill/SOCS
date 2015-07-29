#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 15, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

Simulate the influence of different P0 on the age determination of a cluster
using gyrochronology
'''

import config
import numpy as np
import matplotlib.pyplot as plt

clusterage = 440

def findp0raw(bv, p):
    from gyroage import gyroage
    p0s = np.arange(0.1, 3.4, 0.01)
    ages = []
    for p0 in p0s:
        age = gyroage(bv, p, P0 = p0)
        ages.append(age)
    ages = np.array(ages)
    #plt.plot(p0s, ages)
    #plt.show()
    i = np.argmin(abs(ages - clusterage))
    print i, ages[i], p0s[i]
    return p0s[i]

def findp0(bv, p):
    from scipy.optimize import minimize
    from gyroage import gyroage
    
    def min_func(x, bv, p):
        x = abs(x)
        x = max([x,0.1])
        x = min([x,3.4])
        return abs(gyroage(bv, p, P0=x) - clusterage)
    
    x0 = 1.1
    res = minimize(min_func, x0, method='nelder-mead', args=(bv, p),
                       options={'xtol': 1.0, 'disp': False})
    #print res
    return res['x'][0]
    
def agesim():
    data = np.loadtxt(config.datapath+'periods.txt')
    bv = data[:,1]
    p = data[:,5]
    n = np.shape(p)[0]
    p0 =[]
    for bvi, pi in zip(bv-0.03, p):
        p = findp0raw(bvi, pi)
        p0.append(p)
    
    p0 = np.array(p0)
    print p0
    plt.hist(p0, bins=20)
    plt.title('M48 age = %d' % clusterage)
    plt.xlabel('P$_0$ in days')
    plt.show()
    plt.plot(bv, p0,'o')
    plt.show()

def plot():
    from gyroage import gyroperiod
    data = np.loadtxt(config.datapath+'periods.txt')
    bv = data[:,1]
    p = data[:,5]
    bv400 = np.linspace(0.515, 1.6, num=100)
    P11 = gyroperiod(bv400, clusterage)
    P34 = gyroperiod(bv400, clusterage, P0=3.4)
    P01 = gyroperiod(bv400, clusterage, P0=0.1)
    
    plt.scatter(bv-0.03, p, edgecolor='none', c='k', label='STELLA rotation periods')
    plt.axvline(0.653, color='y', linewidth=3.0)
    
    plt.plot(bv400, P34, 'b--', label='P$_0$=3.4 d')
    plt.plot(bv400, P01, 'g--', label='P$_0$=0.1 d')
    plt.plot(bv400, P11, 'r', linewidth=3.0, label='P$_0$=1.1 d')
    plt.legend()
    
    plt.show()

if __name__ == '__main__':
    #plot()
    agesim()