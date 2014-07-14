'''
Created on Oct 25, 2013

@author: jwe <jweingrill@aip.de>
'''

def filename2starid(filename):
    import os
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0]

def set_good(starid):
    query = """UPDATE m48stars
               SET good=True
               WHERE starid like '%s';""" % starid
    wifsip.execute(query)

def set_bad(starid):
    query = """UPDATE m48stars
               SET good=False
               WHERE starid like '%s';""" % starid
    wifsip.execute(query)


if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    pdfs = glob('/work2/jwe/m48/plots/good/2014*.pdf')
    for p in pdfs:
        #print p
        objid = filename2starid(p)
        print objid,'= good'
        set_good(objid)

    pdfs = glob('/work2/jwe/m48/plots/2014*.pdf')
    for p in pdfs:
        #print p
        starid = filename2starid(p)
        print starid,'= bad'
        set_bad(starid)


    wifsip.commit()
    wifsip.close()
