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

def set_good(objid):
    query = """UPDATE frames
               SET good=True
               WHERE objid like '%s';""" % objid
    wifsip.execute(query)

def set_bad(objid):
    query = """UPDATE frames
               SET good=False
               WHERE objid like '%s';""" % objid
    wifsip.execute(query)


if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    pngs = glob('/work2/jwe/stella/wifsip/m48/rot/science*.png')
    for p in pngs:
        #print p
        objid = filename2objid(p)
        print objid
        set_good(objid)

    pngs = glob('/work2/jwe/stella/wifsip/m48/rot/bad/science*.png')
    for p in pngs:
        #print p
        objid = filename2objid(p)
        print objid
        set_bad(objid)


    wifsip.commit()
    wifsip.close()
