#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 2, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
class ClusterGroup():
    """
    creates a XML group file for a group of frames
    """
    title = 'SOCS'
    abstract = 'Photometric monitoring of open stellar clusters'
    filename = 'group.xml'
    user = 'Weingrill'
    address = 'jweingrill@aip.de'
    moondistance = 30.0
    alttarget = 30.0
    pernight = 25
    timeout = 0
    periodday = 0.015
    stretch = 1.0
    zerofraction = 0.011111
    duration = 1185.00
    
    def __init__(self):
        pass

    def tofile(self, path='./'):
        """
        creates the XML structure and writes it to the given filename
        """
        
        from lxml import etree as ET
        from datetime import datetime as dt
        
        
        def addtext(parent, tag, value):
            """
            adds a simple element of the form e.g. '<User>Weingrill</User>'
            """
            element = ET.SubElement(parent, tag)
            element.text = str(value)
        
        def addconstraint(parent, name, values):
            constraint = ET.SubElement(parent, 'Constraint')
            variable = ET.SubElement(constraint, 'Variable')
            variable.text = name
            for key in values:
                node = ET.SubElement(constraint, key)
                node.text = str(values[key])

        def addconstant(parent, javaclass, name, value):
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
        targetname.text = 'M 67 rot group'

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
        retrymax = ET.SubElement(exception1, 'Retrymax')
        retrymax.text='5'
        addtext(exception1, 'Retry', 'NO_STAR_ON_ACQUIRE')
        addtext(exception1, 'Retry', 'DROP_TARGET')
        exception2 = ET.SubElement(target, 'Exception', {'for': 'this'})
        addtext(exception2, 'Delay', 'NO_DAUGHTER_AVAILABLE')
        
        select = ET.SubElement(target, 'Select')
        addtext(select, 'Requires', 'Roofopen && TelescopeAvailable && WifsipAvailable && InitTelescope')
        
        # - Select

        # calculate the JD for beginning and end of observation
        from astropy.time import Time
        t = Time(now.isoformat(), format='isot', scale='utc')
        print round(t.jd) -0.5
        
        addconstraint(select, 'Jd', {'From': round(t.jd) -0.5, 'To': round(t.jd + 356) -0.5})
        addconstraint(select, 'MoonDistance', {'Min': self.moondistance})
        addconstraint(select, 'AltTarget', {'Min': self.alttarget})
        addconstraint(select, 'SolHeight', {'Max': -16.0})
        
        # -- Merit
        merit = ET.SubElement(select, 'Merit')
        timeslot = ET.SubElement(merit, 'Timeslot', {'class':'stella.merit.PerNightMerit'})
                
        addconstant(timeslot, 'java.lang.String', 'time', 'Time')
        addconstant(timeslot, 'java.lang.Integer', 'pernight', self.pernight)
        addconstant(timeslot, 'java.lang.Boolean', 'endobserve', 'false')
        addconstant(timeslot, 'java.lang.String', 'nightlength', 'NightLength')
        addconstant(timeslot, 'java.lang.Long', 'timeout', self.timeout)
        
        gain = ET.SubElement(merit, 'Gain', {'class':'stella.merit.FixedDelayMerit'})
        
        addconstant(gain, 'java.lang.Double', 'period_day', self.periodday)
        addconstant(gain, 'java.lang.Double', 'stretch', self.stretch)
        addconstant(gain, 'java.lang.Double', 'zerofraction', self.zerofraction)
        # </Merit>
        # </Select>
        addtext(target, 'Duration', self.duration)
        addtext(target, 'Daughter', 'M 67 rot NE')
        setup = ET.SubElement(target, 'Setup', id='priority')
        addtext(setup, 'Instrument', 'SCS')
        addconstant(setup, 'java.lang.Double', 'PriorityCap', 3.0)
        # </Setup>
        objectnode = ET.SubElement(target, 'Object', id='main')
        addtext(objectnode, 'ObjectName', 'HTM NB4AE5F')
        position =  ET.SubElement(objectnode, 'Position')
        #TODO: handle coordinates
        addtext(position, 'RA', 132.825)
        addtext(position, 'Dec', 11.800)
        addtext(position, 'Epoch', 2000.0)
        # </Position>
        # </Object>
        # </Target>
        
        # write the tree to a file
        tree = ET.ElementTree(target)
        doctype = '<!DOCTYPE Target SYSTEM "/stella/home/stella/stella/xml/target.dtd">'
        s =  ET.tostring(tree, encoding="UTF-8", doctype=doctype, pretty_print=True)
        
        f = open(path + self.filename, 'wt')
        f.write(s)
        f.close()
        
if __name__ == '__main__':
    cg = ClusterGroup()
    cg.tofile(path = '/work2/jwe/m67/')
