import numpy as np
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

def read_vpwcs(name,path_data='',base_name='NAME_wcs.txt',hdrid='#'):
    fibid=[]
    ra=[]
    dec=[]
    file=path_data+base_name.replace('NAME',name)
    f=open(file,'r')
    for line in f:
        if not hdrid in line:
            data=line.replace('\n','').split(' ')
            data=list(filter(None,data))
            fibid.extend([int(data[0])-1])
            ra.extend([float(data[1])*3600.0])
            dec.extend([float(data[2])*3600.0])
    f.close()
    fibid=np.array(fibid)
    ra=np.array(ra)
    dec=np.array(dec)
    return fibid,ra,dec

def astrometry_recal(name,path_data='',base_name='NAME_wcs.txt',dra=0,ddec=0,pix_s=0.5,fsc=1.0,returnt=False):
    """
    Recalculate the astrometry for the fibers based on the provided RA and Dec offsets.
    Returns:
        tuple: Updated fiber IDs, RA offsets, and Dec offsets.
    """
    fibid,ra_fib,dec_fib=read_vpwcs(name,path_data=path_data,base_name=base_name)
    ra0t=np.mean(ra_fib-dra)/3600.0
    dec0t=np.mean(dec_fib-ddec)/3600.0
    wt1 = WCS(naxis=2)    
    wt1.wcs.crpix = [1, 1]
    wt1.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
    wt1.wcs.crval = [ra0t,dec0t]
    wt1.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wt1.wcs.radesys = 'ICRS'
    sky_coord = SkyCoord(ra=(ra_fib-dra)/3600.0, dec=(dec_fib-ddec)/3600.0, frame="icrs", unit="deg")           
    x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
    x_pixel=x_pixel*fsc
    y_pixel=y_pixel*fsc
    skycor = pixel_to_skycoord(x_pixel,y_pixel,wt1)
    ran=skycor.ra.value
    decn=skycor.dec.value
    if returnt == False:
        f=open(path_data+base_name.replace('NAME',name+'_New'),'w')
        for i in range(len(fibid)):
            f.write(f'{fibid[i]+1} {ran[i]} {decn[i]}\n')
        f.close()
    else:
        decn=decn*3600.0
        ran=ran*3600.0
        return fibid,ran,decn