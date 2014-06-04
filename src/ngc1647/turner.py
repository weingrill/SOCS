'''
Created on Dec 11, 2013

@author: jwe
'''
class Turner(object):
    def __init__(self):
        """load """
        self.stars = []
    
    def loadtable(self):
        f = open('/work1/jwe/NGC1647/Turner.tsv', 'rt')
        lines = f.readlines()
        f.close()
        
        # skip head
        for l in lines[1:len(lines)-3]:
            sl = l.split('\t')
            self.stars.append({'star': sl[0], 
                               'rem': sl[1], 
                               'V': float(sl[2]),
                               'B-V': float(sl[3]),
                               'U-B': float(sl[4]),
                               'n': int(sl[5]),
                               'code': int(sl[6]),
                               'RA': 0.0,
                               'Dec': 0.0})
    
    def get_coordinates(self):
        import PySimbad
        import astronomy as ast
        for s in self.stars:
            name = s['star']
            print name,
            try:
                coord = PySimbad.SimbadCoord(name)
            except KeyError:
                print 'none'
            else:
                print coord, 
                hms = coord.split(' ')
                s['RA'] = ast.hms2dd(' '.join(hms[1:4]))
                s['Dec'] = ast.dms2dd(' '.join(hms[5:8]))
                print s['RA'], s['Dec']
    
    def create_db(self):
        query = """DROP TABLE IF EXISTS turner;
CREATE TABLE turner( 
 star  varchar(12), 
 rem   varchar(2),
 vmag  real,
 bv    real,
 ub    real,
 n     integer,
 code  integer,
 ra    float,
 dec   float,
 coord point,
 PRIMARY KEY (star,rem));"""
        from datasource import DataSource 
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        wifsip.execute(query)
        wifsip.close()

    def update_db(self):
        from datasource import DataSource 
        import StringIO
        
        wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        columns = ['star','rem','vmag','bv','ub','n','code','ra','dec'] 
        values = ''       
        for s in self.stars:
            values += '%s\t%s\t%.2f\t%.2f\t%.2f\t%d\t%d\t%f\t%f\n' % \
            (s['star'], s['rem'], s['V'], s['B-V'], s['U-B'], s['n'], s['code'],
             s['RA'], s['Dec'])
        f = StringIO.StringIO(values)
        print values
        try:
            wifsip.cursor.copy_from(f,'turner', columns=columns)
        finally:
            wifsip.commit()
        wifsip.close()
        
    
if __name__ == '__main__':
    t = Turner()
    t.loadtable()
    print t.stars[0]
    t.get_coordinates()
    print t.stars[0]
    #t.create_db()
    t.update_db()