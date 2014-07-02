'''
Created on May 3, 2013

@author: jwe <jweingrill@aip.de>

get files from 1647 field 1 in rp
'''
def getframes(targetdir='/work2/jwe/stella/wifsip/m48'):
    from datasource import DataSource
    from subprocess import call

    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    query = """SELECT path, filename
             FROM frames, science
             WHERE (filter LIKE 'V' OR FILTER LIKE 'B')
              AND OBJECT LIKE 'M 48 BVI %%'
              AND frames.objid = science.objid
              order by frames.objid"""
    tab = wifsip.query(query)
    for path, filename in tab:
        print path+'/'+filename
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])

if __name__ == '__main__':
    getframes()