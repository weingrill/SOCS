'''
Created on Oct 25, 2013

@author: jwe <jweingrill@aip.de>
'''

def filename2objid(filename):
    #science20101031A-0009EXP0010.png
    
    i = filename.find('science')
    objid = filename[i+7:]
    objid = objid.replace('EXP','-')
    return objid[:19]

def set_good(objids):
    
    #query = """UPDATE frames
    #           SET good=True
    #           WHERE objid in ('%s');""" % ','.join(objids) 
    query = """UPDATE frames
               SET good=True
               WHERE objid like '%s';""" % objids
    print query
    wifsip.execute(query)

if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    objids = set()
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    pngs = glob('/work2/jwe/stella/wifsip/ngc1647/science*.png')
    for p in pngs:
        #print p
        objid = filename2objid(p)
        print objid
        objids.add(objid)
        set_good(objid)

    wifsip.commit()
    wifsip.close()
