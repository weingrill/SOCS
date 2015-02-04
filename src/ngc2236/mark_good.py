'''
Created on Oct 25, 2013

@author: jwe <jweingrill@aip.de>
'''

def corotid(filename):
    #20131107A-0007-0006#2545(221701040).pdf
    
    i = filename.find('(')
    j = filename.find(')')
    corotid = int(filename[i+1:j])
    return corotid

def set_good(corotid):
    query = """UPDATE ngc2236
               SET good=True
               WHERE corotid = %d;""" % corotid
    wifsip.execute(query)

def set_bad(corotid):
    query = """UPDATE ngc2236
               SET good=False
               WHERE corotid = %d;""" % corotid
    wifsip.execute(query)


if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    pdfs = glob('/work2/jwe/NGC2236/plots/good/*.pdf')
    for p in pdfs:
        cid = corotid(p)
        print p, cid, 'good'
        set_good(cid)
    
    pdfs = glob('/work2/jwe/NGC2236/plots/*.pdf')
    for p in pdfs:
        cid = corotid(p)
        print p, cid, 'bad'
        set_bad(cid)

    wifsip.commit()
    wifsip.close()
