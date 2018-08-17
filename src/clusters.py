'''
Created on Dec 3, 2014

@author: jwe
'''
from datasource import DataSource
from astropy import units as u
from astropy.coordinates import SkyCoord  # @UnresolvedImport
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import _config


class Clusters(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.table = DataSource(database=_config.dbname, user=_config.dbuser, host=_config.dbhost)

    def create_table(self):
        query = """DROP TABLE clusters;
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

        f = open(filename, 'rt')
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
            record = {'name': l[0:18].rstrip(),
                      'ra': c.ra.degree,
                      'dec': c.dec.degree,
                      'class': l[40:43].rstrip(),
                      'diam': nfloat(l[45:52]),
                      'd': nint(l[55:60]),
                      'ebv': nfloat(l[65:70]),
                      'logage': nfloat(l[73:79]),
                      'pmra': nfloat(l[84:90]),
                      'epmra': nfloat(l[92:96]),
                      'pmdec': nfloat(l[100:106]),
                      'epmdec': nfloat(l[108:112]),
                      'nc': nint(l[113:118]),
                      'ref1': l[119:123].rstrip(),
                      'rv': nfloat(l[127:134]),
                      'erv': nfloat(l[138:143]),
                      'n': nint(l[147:150]),
                      'ref2': nint(l[155:159]),
                      'me': nfloat(l[162:168]),
                      'eme': nfloat(l[171:176]),
                      'nme': nint(l[177:180]),
                      'trtyp': l[183:191].rstrip(),
                      'coord': '(%f,%f)' % (c.ra.degree, c.dec.degree)}

            # print record
            def nstr(s):
                if len(str(s).strip()) == 0:
                    return '\\N'
                elif type(s) is str:
                    return str(s)
                else:
                    return str(s)

            valline = '\t'.join([nstr(v) for v in record.values()])
            valline = valline.replace('nan', '\\N')
            print(valline)
            values += valline + '\n'

        columns = record.keys()
        f = StringIO.StringIO(values)
        try:
            cur.copy_from(f, 'clusters', columns=columns)
        finally:
            self.table.commit()

    def query(self):

        query = """SELECT name, ra, dec, ebv, diam from clusters
        WHERE (name like 'NGC %' or name like 'IC %') 
        AND diam>10 AND diam<60
        AND d<1500
        AND dec>-15.0
        AND ebv<0.3
        AND logage>8.0 AND logage<=9.5
        AND abs(rv)>2.0;"""
        result = self.table.query(query)
        names = [r[0] for r in result]
        ra = np.array([r[1] for r in result])
        dec = np.array([r[2] for r in result])
        ebv = np.array([r[3] for r in result])
        diam = np.array([r[4] for r in result])

        mycmap = plt.cm.get_cmap('Reds')
        # mycmap.set_under('w')fig, ax_f = plt.subplots()

        _, ax1 = plt.subplots(figsize=(10, 7))
        plt.scatter(ra / 15., dec, s=diam * 4, cmap=mycmap, c=ebv)
        for rai, deci, iname in zip(ra, dec, names):
            if iname in ['IC 4756', 'NGC 2682', 'NGC 2319', 'NGC 2374', 'NGC 7209', 'NGC 7243', 'NGC 7082', 'NGC 225']:
                horizontalalignment = 'right'
                withdash = None
                dashlength = None
            elif iname in ['NGC 2413']:
                horizontalalignment = 'right'
                withdash = True
                dashlength = 20.0
            else:
                horizontalalignment = 'left'
                withdash = None
                dashlength = None
            if withdash:
                plt.text(rai / 15., deci, iname, withdash=withdash, dashlength=dashlength,
                         horizontalalignment=horizontalalignment)
            else:
                plt.text(rai / 15., deci, iname, horizontalalignment=horizontalalignment)
        # plt.draw()
        ax2 = ax1.twiny()

        print(np.arange(13, 0, -1))
        ax1.set_xlabel('right ascension')
        ax1.set_ylabel('declination')
        ax2.set_xlabel('culmination month')
        ax1.set_xticks(np.arange(24))
        ax1.set_yticks(np.arange(-15, 65, 5))
        ax2.set_xticks(np.arange(0, 13) + 0.4333)
        ax2.set_xticklabels(['9', '8', '7', '6', '5', '4', '3', '2', '1', '12', '11', '10'])
        ax1.set_xlim(24, 0)
        # ax2.set_xlim(12, 0)
        ax1.set_ylim(-15, 65)
        ax1.grid()
        plt.minorticks_on()
        cbar = plt.colorbar()
        # cbar.ax.set_yticklabels(['0','1','2','>3'])
        cbar.set_label('E(B - V)', rotation=270)
        plt.savefig('/work2/jwe/SOCS/plots/cluster_query.pdf')
        plt.close()
        print(len(names))
        print('\n'.join(names))


if __name__ == '__main__':
    cl = Clusters()
    # cl.create_table()
    # cl.fromfile()
    cl.query()
