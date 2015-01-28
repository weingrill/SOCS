'''
Created on Jan 27, 2015

@author: jwe
'''

class M48Phot(object):
    '''
    Import photometry from David James
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource

        self.data=[]
        self.table = DataSource(database='wifsip', user='sro', host='pina.aip.de') 
        
    def fromfile(self, filename):
        f = open(filename)
        lines = f.readlines()
        f.close
        p = filename.find('M48-')
        field = filename[p+4:p+6]
        for l in lines:
            sl = l.split()
            record = {'starid': '%s-%s-%3.3d' % (sl[0],field,int(sl[1])),
                      'ccdx': float(sl[2]),
                      'ccdy': float(sl[3]),
                      'ra': float(sl[10]),
                      'dec': float(sl[11]),
                      'vmag': float(sl[12]),
                      'vmag_err': float(sl[13]),
                      'bv': float(sl[14]),
                      'bv_err': float(sl[15]),
                      'vi': float(sl[16]),
                      'vi_err': float(sl[17])
                      }
            if record['bv']>9: record['bv'] = None 
            if record['bv_err']>9: record['bv_err'] = None 
            if record['vi']>9: record['vi'] = None 
            if record['vi_err']>9: record['vi_err'] = None 
            self.data.append(record)
            print record
    def createtable(self):
        '''
        create table for the photometry
        '''
        if not raw_input('press Y to erase m48phot')=='Y':
            return
        
        query = """
        DROP TABLE IF EXISTS m48phot;
        CREATE TABLE m48phot(
         starid varchar(10),
         ccdx real,
         ccdy real,
         ra double precision,
         dec double precision,
         vmag real default 0,
         vmag_err real,
         bv real,
         bv_err real,
         vi real,
         vi_err real,
         coord point,
         PRIMARY KEY (starid));
        GRANT SELECT ON m48phot TO public;
        CREATE INDEX idx_m48phot_coords ON m48phot USING GIST (circle(coord,0));
        """
        self.table.execute(query)
        print "table 'm48phot' created"
        
    def todatabase(self):
        import StringIO
        
        
        cur = self.table.cursor

        def nstr(s):
            if s is None: return '\N'
            elif type(s) is str: return str(s) 
            else: return str(s)
        
        values = ''

        for record in self.data:
            valline = '\t'.join([nstr(v) for v in record.values()])
            #valline = valline.replace('nan', '\N')
            print valline
            values += valline + '\n'

                
        columns = record.keys()
        f = StringIO.StringIO(values)
        try:
            cur.copy_from(f,'m48phot', columns=columns)
        finally:
            self.table.commit()
        
        
if __name__ == '__main__':
    m48p = M48Phot()
    m48p.createtable()
    for field in range(1,10):
        m48p.fromfile('/work2/jwe/m48/phot/M48-F%d-BVIposRADEC.lis' % field)
    m48p.todatabase()