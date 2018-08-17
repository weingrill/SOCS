#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Mar 1, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun

from astropy.utils import iers
import matplotlib.pyplot as plt

iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL, cache=True)

clustername = 'NGC 6709'
t = '2016-11-01 00:00:00'
ocluster = SkyCoord.from_name(clustername)
izana = EarthLocation(lat=28.301195 * u.deg, lon=-16.509209 * u.deg, height=2096 * u.m)
utcoffset = -1 * u.hour  # Eastern Daylight Time
# time = Time(t) - utcoffset
# m33altaz = ocluster.transform_to(AltAz(obstime=time,location=izana))
# "M33's Altitude = {0.alt:.2}".format(m33altaz)

midnight = Time(t) - utcoffset
delta_midnight = np.linspace(-2, 7, 100) * u.hour
oclusteraltazs = ocluster.transform_to(AltAz(obstime=midnight + delta_midnight, location=izana))

# plt.plot(delta_midnight, oclusteraltazs.secz)
# plt.xlim(-2, 7)
# plt.ylim(1, 4)
# plt.xlabel('Hours from EDT Midnight')
# plt.ylabel('Airmass [Sec(z)]')
# plt.show()

delta_midnight = np.linspace(-12, 12, 1000) * u.hour
times = midnight + delta_midnight
altazframe = AltAz(obstime=times, location=izana)
sunaltazs = get_sun(times).transform_to(altazframe)
oclusteraltazs = ocluster.transform_to(altazframe)

plt.plot(delta_midnight, sunaltazs.alt, color='y', label='Sun')
plt.scatter(delta_midnight, oclusteraltazs.alt, c=oclusteraltazs.az, label=clustername, lw=0, s=8)
# plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
# plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left')
plt.xlim(-12, 12)
plt.xticks(np.arange(13) * 2 - 12)
plt.ylim(0, 90)
plt.xlabel('Hours from EDT Midnight')
plt.ylabel('Altitude [deg]')
plt.title(t)
plt.savefig('/work2/jwe/SOCS/' + clustername + '.pdf')
# plt.show()
plt.close()
