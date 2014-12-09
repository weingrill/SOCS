'''
Created on Dec 3, 2014

@author: jwe
'''

class Clusters(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
    
        self.table = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        
    def create_table(self):
        query="""DROP TABLE clusters;
        CREATE TABLE clusters (
        name   varchar(18),
        ra     double precision,
        dec    double precision,
        class  varchar(3),
        diam   real,
        d      int,
        ebv    real,
        logage real,
        pmra   real,
        epmra  real,
        pmdec  real,
        epmdec real,
        nc     int,
        ref1   varchar(3),
        rv     real,
        erv    real,
        n      int,
        ref2   varchar(4),
        me     real,
        eme    real,
        nme    real,
        trtyp  varchar(8),   
        coord point,
        PRIMARY KEY (name));"""
        self.table.execute(query)
        
    def fromfile(self, filename='/work2/jwe/SOCS/clusters.txt'):
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        import StringIO

        f = open(filename,'rt')
        lines = f.readlines()
        f.close()
        
        def nfloat(s):
            from numpy import nan
            try:
                result = float(s)
            except ValueError:
                return nan
            return result

        def nint(s):
            from numpy import nan
            try:
                result = int(s)
            except ValueError:
                return nan
            return result

        cur = self.table.cursor
        
        values = ''
        
        for l in lines:
            coords = l[18:37] 
            c = SkyCoord(coords, 'icrs', unit=(u.hourangle, u.deg))  # @UndefinedVariable
            record = {'name':   l[0:18].rstrip(),
                      'ra':     c.ra.degree,
                      'dec':    c.dec.degree,
                      'class':  l[40:43].rstrip(),
                      'diam':   nfloat(l[45:52]),
                      'd':      nint(l[55:60]),
                      'ebv':    nfloat(l[65:70]),
                      'logage': nfloat(l[73:79]),
                      'pmra':   nfloat(l[84:90]),
                      'epmra':  nfloat(l[92:96]),
                      'pmdec':  nfloat(l[100:106]),
                      'epmdec': nfloat(l[108:112]),
                      'nc':     nint(l[113:118]),
                      'ref1':   l[119:123].rstrip(),
                      'rv':     nfloat(l[127:134]),
                      'erv':    nfloat(l[138:143]),
                      'n':      nint(l[147:150]),
                      'ref2':   nint(l[155:159]),
                      'me':     nfloat(l[162:168]),
                      'eme':    nfloat(l[171:176]),
                      'nme':    nint(l[177:180]),
                      'trtyp':  l[183:191].rstrip(),   
                      'coord': '(%f,%f)' % (c.ra.degree, c.dec.degree)}
            #print record
            def nstr(s):
                if len(str(s).strip())==0: return '\N'
                elif type(s) is str: return str(s) 
                else: return str(s)
            valline = '\t'.join([nstr(v) for v in record.values()])
            valline = valline.replace('nan', '\N')
            print valline
            values += valline + '\n'

                
        columns = record.keys()
        f = StringIO.StringIO(values)
        try:
            cur.copy_from(f,'clusters', columns=columns)
        finally:
            self.table.commit()

    def query(self):
        import matplotlib.pyplot as plt
        from numpy import array
        
        query = """SELECT name, ra, dec, ebv, diam from clusters
        WHERE (name like 'NGC %' or name like 'IC %') 
        AND diam>10 AND diam<60
        AND d<1000
        AND dec>-10.0
        AND ebv<0.3
        AND logage<=9
        AND abs(rv)>2.0;"""
        result = self.table.query(query)
        names = [r[0] for r in result]
        ra = array([r[1] for r in result])
        dec = array([r[2] for r in result])
        ebv = array([r[3] for r in result])
        diam = array([r[4] for r in result])
        
        plt.scatter(ra/15., dec, s=diam*4, c = ebv)
        for i in zip(ra,dec,names):
            plt.text(i[0]/15., i[1], i[2])
        #plt.draw()  
        plt.xticks()
        plt.xlim(24,0)
        plt.ylim(-10,65)
        plt.show()            
        print len(names)
        print '\n'.join(names)
cl = Clusters()
#cl.create_table()
#cl.fromfile()
cl.query()