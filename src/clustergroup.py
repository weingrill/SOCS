#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 2, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

Global properties of parsed target: {imagetype=object, night=NightRemain, accelerate=1.0, rateretry=10, ratewait=60000, obsbad=false, obsgood=false}

'''
from opencluster import OpenCluster
from astropy.time import Time
from lxml import etree as ET  # @UnresolvedImport
from datetime import datetime as dt
import os
from subprocess import call
import time
        

class ClusterGroup():
    """
    creates a XML group file for a group of frames
    """
    title = 'SOCS'
    abstract = 'Photometric monitoring of open stellar clusters'
    filename = 'group.xml'
    user = 'Weingrill'
    address = 'jweingrill@aip.de'
    #timeout = 0 # should be zero according to tgranzer
    stretch = 1.0
    prioritycap = 3.0
    retrymax = 5
    #zerofraction = 0.011111 # should be low according to tgranzer
    
    daughters = []
    
    
    def __init__(self, opencluster):
        """
        Set basic constraints derived from the OpenCluster class
        """
        if not type(opencluster) is OpenCluster:
            raise TypeError('expecting OpenCluster Object')
        self.ra = opencluster.object['RA']
        self.dec = opencluster.object['Dec']
        self.targetname = opencluster.uname+' group'
        self.filename = self.targetname.lower().replace(' ','_')+'.xml'
        self.objectname = opencluster.objectname
        self.moondistance = opencluster.constraints['MoonDistance.Min']
        self.alttarget = opencluster.constraints['AltTarget.Min']
        self.solheight = opencluster.constraints['SolHeight.Max']
        self._pernight = opencluster.pernight
        self.timeout = opencluster.timeout
        self._periodday = opencluster.mode['period_day']
        self.zerofraction = opencluster.mode['zerofraction']
        self.duration = opencluster.duration
        
        self.startdate = opencluster.startdate
        self.enddate = opencluster.enddate
        assert(self.startdate<self.enddate)
        self.daughters = []
        

    def add_daughter(self, daughtername):
        """
        add daughter fields
        """
        self.daughters.append(daughtername)

    @property
    def pernight(self):
        """
        calculate pernight as global cluster per night times daughterfields
        """
        return self._pernight*len(self.daughters)
    
    @property
    def periodday(self):
        """
        calculate periodday
        """
        if len(self.daughters)>0:
            return self._periodday/len(self.daughters)
        else:
            return self._periodday
    
    def tofile(self, path='./'):
        """
        creates the XML structure and writes it to the given filename
        """
        
        
        def addtext(parent, tag, value):
            """
            adds a simple element of the form e.g. '<User>Weingrill</User>'
            """
            element = ET.SubElement(parent, tag)
            element.text = str(value)
        
        def addconstraint(parent, name, values):
            """
            adds a constraint element to the xml file
            """
            constraint = ET.SubElement(parent, 'Constraint')
            variable = ET.SubElement(constraint, 'Variable')
            variable.text = name
            for key,value in sorted(values.items()):
                node = ET.SubElement(constraint, key)
                node.text = str(value)

        def addconstant(parent, javaclass, name, value):
            """
            adds a constant element to the xml file
            """
            constant = ET.SubElement(parent, 'Constant', {'class':javaclass})
            constantname = ET.SubElement(constant, 'Constantname')
            constantname.text = name
            constantvalue = ET.SubElement(constant, 'Constantvalue')
            constantvalue.text = str(value)

        # prepare submitted timestamp 
        now = dt.utcnow()
        now = now.replace(microsecond = 0) 
        # Target
        target = ET.Element('Target', type='chain', access="enabled", 
                            submitted=now.isoformat()+' UTC')
        
        targetname = ET.SubElement(target, 'TargetName', proposal='cluster.survey')
        targetname.text = self.targetname

        addtext(target, 'Title', self.title)
        addtext(target, 'Abstract', self.abstract)
        addtext(target, 'File', self.filename)
        addtext(target, 'User', self.user)
        
        email = ET.SubElement(target, 'Email')
        addtext(email, 'Address', self.address)
        notify = ET.SubElement(email, 'Notify')
        ET.SubElement(notify, 'Onblock')
        ET.SubElement(notify, 'Onfirstpick')
        ET.SubElement(notify, 'Oncomplete')
        
        addtext(target, 'Institution', 'CORE')
        addtext(target, 'Team', 'Barnes')
        
        ET.SubElement(target, 'History')
        exception1 = ET.SubElement(target, 'Exception', {'for': 'daughter'})
        addtext(exception1, 'Retrymax', self.retrymax)
        addtext(exception1, 'Retry', 'NO_STAR_ON_ACQUIRE')
        addtext(exception1, 'Retry', 'DROP_TARGET')
        
        exception2 = ET.SubElement(target, 'Exception', {'for': 'this'})
        addtext(exception2, 'Delay', 'NO_DAUGHTER_AVAILABLE')
        
        select = ET.SubElement(target, 'Select')
        addtext(select, 'Requires', 'Roofopen && TelescopeAvailable && WifsipAvailable && InitTelescope')
        
        # - Select

        # calculate the JD for beginning and end of observation
        startjd = Time(self.startdate).jd 
        endjd = Time(self.enddate).jd 
        if startjd > endjd:
            raise ValueError('From JD must be earlier than end JD')
        addconstraint(select, 'Jd', {'From': startjd, 'To': endjd})
        if not self.moondistance is None:
            addconstraint(select, 'MoonDistance', {'Min': self.moondistance})
        if not self.alttarget is None:
            addconstraint(select, 'AltTarget', {'Min': self.alttarget})
        addconstraint(select, 'SolHeight', {'Max': self.solheight})
        
        # -- Merit
        merit = ET.SubElement(select, 'Merit')
        timeslot = ET.SubElement(merit, 'Timeslot', {'class':'stella.merit.PerNightMerit'})
                
        addconstant(timeslot, 'java.lang.String', 'time', 'Time')
        addconstant(timeslot, 'java.lang.Integer', 'pernight', self.pernight)
        addconstant(timeslot, 'java.lang.Boolean', 'endobserve', 'false')
        addconstant(timeslot, 'java.lang.String', 'nightlength', 'NightLength')
        addconstant(timeslot, 'java.lang.Long', 'timeout', self.timeout)
        
        gain = ET.SubElement(merit, 'Gain', {'class':'stella.merit.FixedDelayMerit'})
        
        if not self.periodday>0.0:
            raise ValueError('period_day must be greater than zero')
        
        addconstant(gain, 'java.lang.Double', 'period_day', self.periodday)
        addconstant(gain, 'java.lang.Double', 'stretch', self.stretch)
        
        if not self.zerofraction>0.0:
            raise ValueError('zerofraction must be greater than zero')
        
        addconstant(gain, 'java.lang.Double', 'zerofraction', self.zerofraction)
        # </Merit>
        # </Select>
        addtext(target, 'Duration', self.duration)
        
        for daughter in self.daughters:
            addtext(target, 'Daughter', daughter )
        setup = ET.SubElement(target, 'Setup', id='priority')
        addtext(setup, 'Instrument', 'SCS')
        addconstant(setup, 'java.lang.Double', 'PriorityCap', self.prioritycap)
        # </Setup>
        objectnode = ET.SubElement(target, 'Object', id='main')
        addtext(objectnode, 'ObjectName', self.objectname)
        position =  ET.SubElement(objectnode, 'Position')
        if self.ra<0.0 or self.ra>360.0:
            raise ValueError('R.A. coordinates out of bounds')
        addtext(position, 'RA', self.ra)
        if self.dec<-90.0 or self.dec>90.0:
            raise ValueError('Dec. coordinates out of bounds')
        addtext(position, 'Dec', self.dec)
        addtext(position, 'Epoch', 2000.0)
        addtext(position, 'Equinox', 2000.0)
        addtext(position, 'V', 22.0)
        addtext(position, 'B-V', 0.0)
        
        # </Position>
        # </Object>
        # </Target>
        
        # write the tree to a file
        tree = ET.ElementTree(target)
        doctype = '<!DOCTYPE Target SYSTEM "/stella/home/stella/stella/xml/target.dtd">'
        s =  ET.tostring(tree, encoding="UTF-8", doctype=doctype, pretty_print=True)
        # remember the path for later use of transfer method
        self.filepath = path
        with open(os.path.join(path,self.filename), 'wb') as f:
            f.write(s)
        
        
    def transfer(self, path=None):
        """
        a plain copy from Opencluster.transfer method
        """
        
        # if path has not been set use the one, where the file has been created. 
        if path is None:
            path = self.filepath
        source = os.path.join(path, self.filename)
        #xml goes dirctly to submit
        target='sro@habanero:/stella/home/www/uploads/weingrill/save/'
        time.sleep(1) # otherwise submit.jnlp gets confused
        print('scp %s %s' % (source, target))
        call(['/usr/bin/scp', source, target])
        print(os.path.dirname(source))
        _, filename = os.path.split(source)
        call(['/usr/bin/ssh', 'operator@ciruelo', 'bin/groupsubmit.sh %s' % filename])
        
