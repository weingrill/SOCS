'''
Created on Apr 12, 2013

@author: jwe <jweingrill@aip.de>
'''

def do_uvby(transfer=False, plot=False):
    from opencluster import OpenCluster
    
    m67 = OpenCluster(objectname='M 67', 
                          uname='M 67 uvby', 
                          obsmode='uvby')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'
    m67.tofile('/work2/jwe/m67')
    if transfer:
        m67.transfer()

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(1.)
        m67.tycho()
        m67.plot()
        plt.show()
    
def do_hahb(transfer=False):
    from opencluster import OpenCluster
    
    m67 = OpenCluster(objectname='M 67', 
                          uname='M 67 hahb', 
                          obsmode='hahb')
    m67.title = 'SOCS'
    m67.abstract = 'Photometric monitoring of open stellar clusters'
    m67.tofile('/work2/jwe/m67')
    if transfer:
        m67.transfer()

        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    do_uvby(transfer=True)
    do_hahb(transfer=True)
