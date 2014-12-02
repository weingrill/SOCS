'''
Created on Sep 16, 2014

@author: jwe

global configurationfile for NGC2236
'''
import os

projectpath = '/work2/jwe/NGC2236/'

resultpath = projectpath+'results/'
plotpath = projectpath+'plots/'
datapath = projectpath+'data/'
lightcurvespath = projectpath+'lightcurves/'

if not os.path.exists(projectpath): os.mkdir(projectpath)
if not os.path.exists(resultpath): os.mkdir(resultpath)
if not os.path.exists(plotpath): os.mkdir(plotpath)
if not os.path.exists(datapath): os.mkdir(datapath)
if not os.path.exists(lightcurvespath): os.mkdir(lightcurvespath)