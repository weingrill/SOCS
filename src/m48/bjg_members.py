'''
Created on Aug 7, 2014

@author: jwe
'''
from astronomy import hms2dd,dms2dd
from datasource import DataSource

wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')

def setmember(ra, dec, membership):
    if membership:
        mem='TRUE'
    else: mem='FALSE'
    query = '''UPDATE m48stars
    SET member=%s
    WHERE circle(point(%f,%f),0) <@ circle(coord,3./3600)''' % (mem, ra, dec)
    wifsip.execute(query)
    

f = open('/work2/jwe/m48/data/table3.dat')
lines = f.readlines()
f.close

for l in lines:
    sl = l.split()
    num = sl[0]
    rah,ram,ras = sl[1:4]
    ra = hms2dd([rah,ram,ras])
    ded,dem,des = sl[4:7]
    dec = dms2dd([ded,dem,des])
    if sl[19]=='M': member = True
    elif sl[19]=='NM': member = False
    print num, ra,dec, member
    setmember(ra, dec, member)