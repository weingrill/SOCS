#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 11, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
import pylab as plt
from datasource import DataSource

wifsip = DataSource(database='wifsip', user='sro', host='oldpina.aip.de')
query = """SELECT bv, period 
            FROM m48stars 
            WHERE good;"""

data = wifsip.query(query)

bv = np.array([d[0] for d in data])
period = np.array([d[1] for d in data])

from gyroage import gyroperiod

bv400 = np.linspace(0.44, 1.6, num=100)
P11 = gyroperiod(bv400, 400)
P34 = gyroperiod(bv400, 400, P0=3.4)
P01 = gyroperiod(bv400, 400, P0=0.1)
P800 = gyroperiod(bv400, 800, P0=1.1)

plt.scatter(bv-0.03, period, edgecolor='none', c='k', label='STELLA rotation periods')
plt.axvline(0.653, color='y', linewidth=3.0)

plt.plot(bv400, P34, 'b--', label='400 Myr model P$_0$=3.4 d')
#plt.plot(bv400, P01, 'g--', label='400 Myr model P$_0$=0.1 d')
plt.plot(bv400, P11, 'r', linewidth=3.0, label='400 Myr model P$_0$=1.1 d')
plt.plot(bv400, P800, 'g', label='800 Myr model P$_0$=1.1 d')
plt.legend()


plt.legend(loc='lower center', fontsize='small')
 
plt.xlabel('stellar color (B - V)$_0$')
plt.ylabel('stellar rotation period [days]')
plt.ylim(0.0, 14.0)
plt.xlim(0.4, 1.6)
plt.minorticks_on()
plt.grid()
plt.savefig('/home/jwe/Documents/Talks/Heraeus2015/m48cpd.svg',transparent=True)
plt.savefig('/home/jwe/Documents/Talks/Heraeus2015/m48cpd.png',transparent=True)
plt.close()
