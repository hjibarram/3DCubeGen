import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

def median_a(x,lw=5,lower=10000,wave=[]):
    if len(wave) > 0:
        index=np.where(wave < lower)[0]
        index2=np.where(wave >= lower)[0]
        x1=np.copy(x)
        x=x[index]
    x_n=np.zeros(len(x))
    for i in range(0, len(x)):
        if i <= lw:
            x_d=x[0:lw]
            #x_d=reject_outliers(x_d)
            x_n[i]=np.nanmean(x_d)
        if i >= len(x)-lw:
            x_d=x[len(x)-lw:len(x)]
            #x_d=reject_outliers(x_d)
            x_n[i]=np.nanmean(x_d)
        if i > lw and i < len(x)-lw:
            x_d=x[i-lw:i+lw]
            #x_d=reject_outliers(x_d)
            x_n[i]=np.nanmean(x_d) 
    if len(wave) > 0:
        x1[index]=x_n
        x1[index2]=x_n[-1]
        x_n=x1
    return x_n


def read_explist(fname='Orion'):
    ft=open(fname,'r')
    tileid=[]
    tilegp=[]
    mjd=[]
    expn=[]
    for line in ft:
        if not '#' in line: 
            if not 'TILE' in line:
                data=line.replace('\n','').replace(' ','').split(',')
                data=list(filter(None,data))
                tileid.extend([data[0]])
                if data[0] == '11111':
                    tg='0011XX'
                else:
                    tg=data[0][0:4]+'XX'
                tilegp.extend([tg])
                mjd.extend([data[1]])
                expn.extend([data[2]])#{:0>8}".format(expT[i]
    return expn,mjd,tileid,tilegp


def read_lvlist(fname='LV_obs.csv'):
    ft=open(fname,'r')
    name=[]
    dit=[]
    for line in ft:
        if not 'target' in line: 
            if not '#' in line:
                data=line.replace('\n','').replace(' ','').split(',')
                data=list(filter(None,data))
                name.extend([data[0]])
                dit.extend([int(data[3])])
    return name,dit

def read_spall(object_name='Orion',redux_dir='',version='1.0.2.dev0'):
    from pandas import HDFStore
    try:
        hdu_list = fits.open('LVM.fits')
    except:
        print('copy and save as fits the dither summary of https://lvm-viki.lco.cl/progress.html to create the LVM.fits file') 
        return
    table_hdu = hdu_list[1]
    table_data = table_hdu.data
    target=table_data.field('target')
    tile_id=table_data.field('tile_id')
    jd=table_data.field('jd')
    mjd=np.copy(jd)
    for i in range(0, len(jd)):
        mjd[i]=jd[i]-2400000
    store0 =HDFStore(redux_dir+'/'+version+'/drpall-'+version+'.h5')
    #store0['summary'].to_csv('outputFileForTable2.csv')
    dset=store0['summary']
    mjd_t=np.array(dset['mjd'])
    tile_g=np.array(dset['tilegrp'])
    tile_idt=np.array(dset['tileid'])
    exp_n=np.array(dset['expnum'])
    store0.close()
    idt=np.zeros(len(exp_n))
    for i in range(0, len(tile_idt)):
        for j in range(0, len(tile_id)):
            if tile_idt[i] == tile_id[j]:
                if target[j] == object_name:
                    idt[i]=1
    nt=np.where(idt == 1)
    expF=list()
    expT=exp_n[nt]
    mjdT=mjd_t[nt]
    tileT=tile_idt[nt]
    tilegF=list(tile_g[nt])
    expF=[]
    mjdF=[]
    tileF=[]
    for i in range(0, len(expT)):
        expF.extend(["{:0>8}".format(expT[i])])
        mjdF.extend([str(mjdT[i])])
        tileF.extend([str(tileT[i])])
    return expF,mjdF,tileF,tilegF
            
def get_apertures(file):
    ra=[]
    dec=[]
    rad=[]
    colr=[]
    namet=[]
    l1=[]
    l2=[]
    th=[]
    typ=[]
    f=open(file,'r')
    ct=1
    for line in f:
        if not 'Region' in line and not 'fk5' in line and not 'global' in line:
            if 'circle' in line:
                data=line.replace('\n','').replace('circle(','').replace('") # color=',' , ').replace(' width=',' , ').replace(' text={',' , ').replace('}',' ')
                data=data.split(',')
                data=list(filter(None,data))
                #print(data)
                ra.extend([data[0]])
                dec.extend([data[1]])
                rad.extend([float(data[2])])
                colr.extend([data[3].replace(' ','')])
                try:
                    namet.extend([data[5].replace(' ','')])
                except:
                    namet.extend([str(int(ct))])
                l1.extend([np.nan])
                l2.extend([np.nan])
                th.extend([np.nan])    
                typ.extend(['circle'])
            if 'box' in line:
                data=line.replace('\n','').replace('box(','').replace(') # color=',' , ').replace(' width=',' , ').replace(' text={',' , ').replace('}',' ')
                data=data.split(',')
                data=list(filter(None,data))
                ra.extend([data[0]])
                dec.extend([data[1]])
                l1.extend([float(data[2].replace('"',''))])
                l2.extend([float(data[3].replace('"',''))])
                th.extend([float(data[4])])
                colr.extend([data[5].replace(' ','')])
                try:
                    namet.extend([data[7].replace(' ','')])
                except:
                    namet.extend([str(int(ct))])
                rad.extend([np.nan])    
                typ.extend(['box'])
            ct=ct+1
    ra=np.array(ra)
    dec=np.array(dec)
    rad=np.array(rad)
    colr=np.array(colr)
    namet=np.array(namet)
    typ=np.array(typ)
    l1=np.array(l1)
    l2=np.array(l2)
    th=np.array(th)
    return ra,dec,rad,l1,l2,th,colr,namet,typ

def extract_spec(spec,hdr,ra='',dec='',rad=1.5,pix=0.35,avgra=False):
    sky1=SkyCoord(ra+' '+dec,frame=FK5, unit=(u.hourangle,u.deg))
    val1=sky1.ra.deg
    val2=sky1.dec.deg
    wcs = WCS(hdr)
    wcs=wcs.celestial
    ypos,xpos=skycoord_to_pixel(sky1,wcs)
    print(xpos,ypos,'POS Pixel')
    val1=sky1.to_string('hmsdms')
    print(val1,'RA1,DEC1')
        
    nz,nx,ny=spec.shape
    radis=np.zeros([nx,ny])
    for i in range(0, nx):
        for j in range(0, ny):
            x_n=i-xpos
            y_n=j-ypos
            r_n=np.sqrt((y_n)**2.0+(x_n)**2.0)*pix
            radis[i,j]=r_n
    single_T=np.zeros(nz)
    nt=np.where(radis <= rad)
    for i in range(0, nz):
        tmp=spec[i,:,:]
        tmp[np.where(tmp <= 0)]=np.nan
        if avgra:
            single_T[i]=np.nanmean(tmp[nt])
        else:
            single_T[i]=np.nansum(tmp[nt])
        
    crpix=hdr["CRPIX3"]
    try:
        cdelt=hdr["CD3_3"]
    except:
        cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    wave_f=(crval+cdelt*(np.arange(nz)+1-crpix))*1e10
    
    return wave_f,single_T,xpos,ypos

def sycall(comand):
    import os
    linp=comand
    os.system(comand)

def rotate(xx,yy,angle):
    # rotate x and y cartesian coordinates by angle (in degrees)
    # about the point (0,0)
    theta = -1.*angle * np.pi / 180. # in radians
    xx1 = np.cos(theta) * xx - np.sin(theta) * yy
    yy1 = np.sin(theta) * xx + np.cos(theta) * yy
    return xx1, yy1

def make_radec(xx0,yy0,ra,dec,pa):
    xx, yy = rotate(xx0,yy0,pa)
    ra_fib = ra*3600.0 + xx/np.cos(dec*np.pi/180.) 
    dec_fib = dec*3600.0 - yy 
    return ra_fib, dec_fib    