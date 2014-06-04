'''
Created on Mar 18, 2014

@author: jwe

import n2281.pos into the database
'''

class Ngc2281Ref(object):
    '''
    classdocs
    '''


    def __init__(self, filename):
        '''
        Constructor
        '''
        self.filename = filename
        
    def doimport(self):
        '''import the file and shove it into the database'''
        from datasource import DataSource
        from astronomy import hms2dd, dms2dd
        from StringIO import StringIO
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        
        # skip header information and compendium
        lines = lines[71:247]
        print lines[0]
        values = ''
        for l in lines:
            id = l[0:9].strip()
            ra = hms2dd(l[9:19].strip())
            dec = dms2dd(l[20:31].strip())
            s = l[32:33]
            try:
                ucac2 = int(l[34:42])
            except ValueError:
                ucac2 = None
            tyc2 = l[43:54]
            v = float(l[56:61])
            remarks = l[62:].rstrip()
            print id,ra,dec,s,ucac2,tyc2,v,remarks
            values += '%s\t%f\t%f\t%s\t%s\t%s\t%f\t%s\n' % (id,ra,dec,s,ucac2,tyc2,v,remarks)
            values = values.replace('None', '\N') 
        f = StringIO(values)
        print values
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        wifsip.cursor.copy_from(f, 'ngc2281ref', 
                                columns=('id','ra','dec','s','ucac2','tyc2','vmag', 'remarks'))
        wifsip.commit()
        wifsip.close()

    def update(self, filename):
        '''update the V and B-V according to Pesch 1961'''
        from datasource import DataSource
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        lines = lines[1:]
        wifsip = DataSource(host = 'pina', database = 'wifsip', user = 'sro')
        for l in lines:
            ls = l.split()
            star = ls[2]
            v = ls[3]
            bv = ls[4]
            ub = ls[5]
            n = ls[6]
            code = ls[7]
            print star,v,bv,ub,n,code
            query = "UPDATE ngc2281ref SET vmag=%s, bv=%s WHERE id like '%s'" % (v,bv,star)
            wifsip.execute(query)
        wifsip.close() 
        
if __name__ == '__main__':
    ref = Ngc2281Ref('/work2/jwe/NGC2281/n2281.pos')
    #ref.doimport()
    #ref.update('/work2/jwe/NGC2281/pesch1961.txt')