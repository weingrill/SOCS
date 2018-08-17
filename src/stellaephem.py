'''
Created on Apr 10, 2013

@author: jwe <jweingrill@aip.de>
'''

import matplotlib.pyplot as plt

import ephem
import datetime
import numpy as np
from astronomy import airmass

stella = ephem.Observer()
stella.lon, stella.lat = '-16.50925', '28.301215'  # Izana
# stella.lon, stella.lat = '13.104659', '52.404963' #AIP
sun, moon = ephem.Sun(), ephem.Moon()  # @UndefinedVariable

stella.pressure = 0
stella.horizon = '-0:34'
stella.elevation = 2000

ic4756 = ephem.readdb('IC 4756,f|O,18:39: 0,+05:27, 5.,2000,3120')
ic4756.compute()

print('Moonrise:', stella.previous_rising(moon))
print('Moonset: ', stella.next_setting(moon))
print('Sunrise: ', stella.previous_rising(sun))
print('Sunset:  ', stella.next_setting(sun))

today = datetime.datetime.today()
dt = datetime.timedelta(days=14)
today += dt

sun_alt = np.empty(24)
moon_alt = np.empty(24)
minutes = np.arange(24 * 60)
ic4756_alt = np.empty(24)
for m in minutes:
    stella.date = today
    sun.compute(stella)
    moon.compute(stella)
    ic4756.compute(stella)
    sun_alt[h] = float(sun.alt)
    moon_alt[h] = float(moon.alt)
    ic4756_alt[h] = float(ic4756.alt)

fig = plt.figure()
ax_h = fig.add_subplot(111)

ax_h.set_ylim(0, 90)
ax_h.set_xlim(0, 24)
# code was broken, inserted the following line:
hours = np.arange(24)
ax_airmass = ax_h.twinx()

ax_h.set_xticks(hours)
heights = ax_h.get_yticks()
am = airmass(heights)
aml = ['%.2f ' % a for a in am]
ax_airmass.set_ylim(0., 90.)
ax_airmass.set_yticklabels(aml)
ax_h.grid()
ax_h.plot(hours, sun_alt * 180.0 / np.pi, 'yo')
ax_h.plot(hours, moon_alt * 180.0 / np.pi, 'go')
ax_h.plot(hours, ic4756_alt * 180.0 / np.pi, 'k')

ax_h.set_xlabel("hours")
ax_h.set_ylabel("height (degrees)")
ax_airmass.set_ylabel("airmass")

plt.draw()
plt.show()
