#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 14, 2014

@author: jwe
'''

if __name__ == '__main__':
    from functions import scaleto
    from PIL import Image  # @UnresolvedImport
    import numpy as np
    import pyfits
    import pylab as plt
    
    fitsfile = '/work2/jwe/stella/wifsip/m48/coadd.fits'
    fileout = '/work1/jwe/Dropbox/Documents/Poster/AIP2014/m48.jpg'
    
    try:
        hdulist = pyfits.open(fitsfile)
        hdr = hdulist[0].header       
        img = hdulist[0].data
    except IOError:
        exit()
    finally:
        hdulist.close()
        
    fimg = np.flipud(img)
    background = np.mean(fimg)
    print background
    cimg = np.clip(fimg, background, 65535)
    limg = np.log10(cimg)
    simg = scaleto(limg,[255.0,0.0])
    print
    plt.hist(simg.flatten(), bins=256, log=True)
    plt.show()
    a = np.rint(simg)
    b = a.astype('uint8')
    im = Image.fromarray(b)
    im.save(fileout)                    
