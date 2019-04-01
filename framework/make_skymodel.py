import pyfits as pf
import numpy as np
import sys

base_fitsimage = sys.argv[1]
freq = float(sys.argv[2]) # Hz
pixscale = float(sys.argv[3]) # arcsec

def oneGaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, bpa, offset):
    sigma_y = sigma_x    
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(bpa)**2)/(2*sigma_x**2) + (np.sin(bpa)**2)/(2*sigma_y**2)
    b = -(np.sin(2*bpa))/(4*sigma_x**2) + (np.sin(2*bpa))/(4*sigma_y**2)
    c = (np.sin(bpa)**2)/(2*sigma_x**2) + (np.cos(bpa)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( -1.0* (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g #.ravel()



d = np.squeeze(pf.getdata(base_fitsimage))
hdr = pf.getheader(base_fitsimage)
hdr['CDELT1'] = -pixscale/3600.
hdr['CDELT2'] = pixscale/3600.
hdr['CRVAL4'] = freq 
d *= 0
subsize = int((len(d) )/2.)
y=range(len(d[0]))
x=range(len(d[0]))
x, y = np.mgrid[-subsize:subsize, -subsize:subsize]

#oneGaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, bpa, offset)
#g1 = oneGaussian((x,y),1, 0,0,5,5,0,0)
#g2 = oneGaussian((x,y),1, 0,0,20,20,0,0)
#d = (g1/g1.sum()) + (g2/g2.sum()) # two int flux 1 Jy srcs

d[514,514] = 1 # 1 Jy point source

pf.writeto(base_fitsimage.replace('.fits','_ptsrc_gauss.fits'),d.reshape(1,1,len(d),len(d)),header=hdr,clobber=True)
