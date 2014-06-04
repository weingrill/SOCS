'''
Created on May 3, 2013

@author: jwe <jweingrill@aip.de>

get files from 1647 field 1 in rp
'''
def getframes(targetdir='/work2/jwe/stella/wifsip/ngc1647'):
    from datasource import DataSource
    from subprocess import call

    
    wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
    query = """select path, filename
             from frames, science
             where filter like 'rp'
              and object like 'NGC 1647 field 1'
              and expt = 62.5 and frames.objid = science.objid"""
    tab = wifsip.query(query)
    for path, filename in tab:
        print path+'/'+filename
        call(['scp', 'sro@pina.aip.de:'+path+'/'+filename+'.fitz ',targetdir])

if __name__ == '__main__':
    getframes()