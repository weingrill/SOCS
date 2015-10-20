'''
Created on Oct 25, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
def filename2starid(filename):
    import os
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0].rstrip('(0)')

def set_good(starid):
    query = """UPDATE ngc6633
               SET good=True
               WHERE starid like '%s';""" % starid
    wifsip.execute(query)

if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
    pdfs = glob(config.plotpath+'good/*.pdf')
    for p in pdfs:
        #print p
        objid = filename2starid(p)
        print objid,'= good'
        set_good(objid)

    wifsip.commit()
    wifsip.close()
