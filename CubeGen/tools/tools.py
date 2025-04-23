import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d
from astropy.convolution import convolve,Gaussian2DKernel
import CubeGen
import CubeGen.tools.kernel as kernel 
import os
import yaml
from multiprocessing.pool import ThreadPool
from scipy.spatial.distance import pdist
from tqdm.notebook import tqdm

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


def read_explist(fname='Orion',path=''):
    ft=open(path+fname,'r')
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

def simpson_r(f,x,i1,i2,typ=0):
    n=(i2-i1)*1.0
    if n % 2:
        n=n+1.0
        i2=i2+1
    b=x[i2]
    a=x[i1]
    h=(b-a)/n
    s= f[i1]+f[i2]
    n=int(n)
    dx=b-a
    for i in range(1, n, 2):
        s += 4 * f[i1+i]
    for i in range(2, n-1, 2):
        s += 2 * f[i1+i]
    if typ == 0:
        return s*h/3.0
    if typ == 1:
        return s*h/3.0/dx

def get_narrwband(wave, lo=6563,dw=10.0,sig=1.0,alpha=1): 
    y1=(1+np.exp(-((wave-(lo-dw/2))/sig)))**(-alpha)
    y2=(1+np.exp(((wave-(lo+dw/2))/sig)))**(-alpha) 
    y=y1*y2
    return y

def band_spectra(wave_s,pdl_flux,k=5,zt=0):
    dir=os.path.join(CubeGen.__path__[0], 'legacy')+'/'
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    wave_s=wave_s/(1+zt)
    nt=np.where((np.isfinite(pdl_flux) == True) & (pdl_flux > 0))
    pdl_flux=pdl_flux[nt]
    wave_s=wave_s[nt]
    [nw]=pdl_flux.shape
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    file=['LICK_LICK.3734_63.dat','LaSilla_DFOSCOIII_filter.dat','LICK_LICK.6563_100.dat','LICK_LICK.6710_100.dat','SLOAN_SDSS.g.dat','SLOAN_SDSS.r.dat','SLOAN_SDSS.i.dat','Palomar_POSS.Red.dat']
    band=['OII','OIII','HI','SII','g','r','i','R']
    zerop=[3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    f=open(dir+file[k],'r')
    wave=[]
    trans=[]
    for line in f:
        data=line.replace('/n','').split(' ')
        data=list(filter(None,data))
        if len(data) > 1:
            wave.extend([float(data[0])])
            trans.extend([float(data[1])])
    f.close()
    d_wave=np.zeros(len(wave))
    for kk in range(1,len(wave)):
        d_wave[kk]=wave[kk]-wave[kk-1]
    d_wave[0]=d_wave[1]
    trans=np.array(trans)
    wave=np.array(wave)
    if 'OII' == band[k]:
        trans=get_narrwband(wave, lo=3728.48,dw=8)#doublet
    if 'OIII' == band[k]:
        trans=get_narrwband(wave, lo=5008.24,dw=5) 
    if 'HI' == band[k]:
        trans=get_narrwband(wave, lo=6564.61,dw=5) 
    if 'SII' == band[k]:
        trans=get_narrwband(wave, lo=6732.67,dw=5)    
    spec=np.copy(pdl_flux)
    if np.nansum(spec) != 0:
        spec1=interp1d(wave_s, spec,kind='linear',bounds_error=False,fill_value=0.)(wave)
        flux_t=spec1*trans*d_wave
        f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
        f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2,typ=0)#/simpson_r(trans,wave,0,len(wave)-2,typ=1)
        if f_fin <= 0:
            photo_a=-2.5*np.log10(1e-10)#14
        else:
            photo_a=-2.5*np.log10(f_fin)
        photo_c=f_fin*zerop[k]*jans
        photo_b=f_fi2
    else:
        photo_a,photo_b,photo_c=[0,0,0]
    return photo_a,photo_b,photo_c,band[k]

def twoD_interpolB(x,y,x1,x2,x3,y1,y2,y3,z1,z2,z3):
    a=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
    b=(z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
    c=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
    d=-(a*x1+b*y1+c*z1)
    z=-(a*x+b*y+d)/c
    if np.isscalar(z) == False:
        #nt=np.where(np.isfinite(z) == False)
        nt = np.where(~np.isfinite(z))
        if nt[0].size > 0:
            z[nt]=0.0
    else:
        if np.isfinite(z) == False:
            z=0.0
    return z        
    
def map_interpolB(cube,x,y,nxt=10,nyt=10):
    xpos0=np.int(np.round(x))
    ypos0=np.int(np.round(y))
    map=cube[xpos0-nxt:xpos0+nxt,ypos0-nyt:ypos0+nyt]
    nx,ny=map.shape
    rp=np.zeros([nx,ny])
    x_t=np.arange(nx)
    y_t=np.arange(ny)
    rp[:,:]=np.nan
    for i in range(0, nx):
        for j in range(0, ny):
            rp[i,j]=np.sqrt((i+xpos0-nxt-x)**2.0+(j+ypos0-nyt-y)**2.0)
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x1=min_in[0]
    y1=min_in[1]
    rp[x1,y1]=1000.0
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x2=min_in[0]
    y2=min_in[1]
    rp[x2,y2]=1000.0
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x3=min_in[0]
    y3=min_in[1]
    rp[x3,y3]=1000.0
    z1=map[x1,y1]
    z2=map[x2,y2]
    z3=map[x3,y3]
    x1=x1+xpos0-nxt
    y1=y1+ypos0-nyt
    x2=x2+xpos0-nxt
    y2=y2+ypos0-nyt
    x3=x3+xpos0-nxt
    y3=y3+ypos0-nyt
    z=twoD_interpolB(x,y,x1,x2,x3,y1,y2,y3,z1,z2,z3)
    return z   

def cube_interpolB(cube,x,y):
    nz,nx,ny=cube.shape
    rp=np.zeros([nx,ny])
    x_t=np.arange(nx)
    y_t=np.arange(ny)
    rp[:,:]=np.nan
    for i in range(0, nx):
        for j in range(0, ny):
            rp[i,j]=np.sqrt((i-x)**2.0+(j-y)**2.0)
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x1=min_in[0]
    y1=min_in[1]
    rp[x1,y1]=1000.0
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x2=min_in[0]
    y2=min_in[1]
    rp[x2,y2]=1000.0
    min_in=np.unravel_index(np.nanargmin(rp), (nx,ny))
    x3=min_in[0]
    y3=min_in[1]
    rp[x3,y3]=1000.0
    z1=cube[:,x1,y1]
    z2=cube[:,x2,y2]
    z3=cube[:,x3,y3]
    z=twoD_interpolB(x,y,x1,x2,x3,y1,y2,y3,z1,z2,z3)
    return z        
    

def cube_interpol(cube,x,y):
    nz,nx,ny=cube.shape
    val_out=np.zeros(nz)
    for i in range(0, nz):
        #print(i,nz)
        mapt=np.copy(cube[i,:,:])
        zt=twoD_interpol(mapt,x,y)
        val_out[i]=zt
    return val_out    

def extract_segm(hdr,l1=12,l2=12,ra='',dec='',dx=0):
    sky1=SkyCoord(ra+' '+dec,frame=FK5, unit=(u.hourangle,u.deg))
    ra_deg=sky1.ra.deg
    dec_deg=sky1.dec.deg
    ra1=ra_deg-l1/2./3600.
    ra2=ra_deg+l1/2./3600.
    dec1=dec_deg-l2/2./3600.
    dec2=dec_deg+l2/2./3600.
    sky00=SkyCoord(ra1,dec1,frame=FK5, unit=(u.deg,u.deg))
    sky11=SkyCoord(ra2,dec2,frame=FK5, unit=(u.deg,u.deg))
    wcs = WCS(hdr)
    wcs=wcs.celestial
    ypos,xpos=skycoord_to_pixel(sky1,wcs)
    ypos00,xpos00=skycoord_to_pixel(sky00,wcs)
    #print(ypos,xpos,sky1) 
    #print(ypos00,xpos00,sky00)
    xpos00=np.int(np.round(xpos00))
    ypos00=np.int(np.round(ypos00))
    ypos11,xpos11=skycoord_to_pixel(sky11,wcs)
    xpos11=np.int(np.round(xpos11))
    ypos11=np.int(np.round(ypos11))    
    
    return xpos00,xpos11,ypos00-dx,ypos11-dx

def interpolate_matrix(matrix_input,nt=4,ne=2,verbose=False,smoth=True):
    nx,ny=matrix_input.shape
    nx1=int(nx*nt)
    ny1=int(ny*nt)
    matrix_new=np.zeros([nx1,ny1])
    if verbose:
        pbar=tqdm(total=ny1)
    dxt=(nx-ne*2)/float(nx1)
    dyt=(ny-ne*2)/float(ny1)
    xpos=ne
    for i in range(0, nx1):
        xpos=xpos+dxt
        ypos=ne
        for j in range(0, ny1):
            ypos=ypos+dyt
            val=map_interpolB(matrix_input,xpos,ypos,nxt=ne,nyt=ne)
            matrix_new[i,j]=val 
        if verbose:
            pbar.update(1)
    if verbose:
        pbar.close()   
    if smoth:    
        PSF=Gaussian2DKernel(x_stddev=nt/2,y_stddev=nt/2)#4
        matrix_new=convolve(matrix_new, PSF)
    return matrix_new

def read_config_file(file):
    try:
        with open(file, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return data
    except:
        print('Config File not found')
        return None

def kernel_pipe(type='b'):
    #16-a,8-b,4-c,2-d,1-e,1/2-f,1/4-g,1/8-h,1/16-i,1/32-i
    if type == 'a':
        valt=16
    if type == 'b':
        valt=8
    if type == 'c':
        valt=4
    if type == 'd':
        valt=2
    if type == 'e':
        valt=1
    if type == 'f':
        valt=1/2
    if type == 'g':
        valt=1/4
    if type == 'h':
        valt=1/8
    if type == 'i':
        valt=1/16
    if type == 'j':
        valt=1/32
    sigm_s=17.6/2*32*valt
    pix_s=17.6/2*0.75*32*valt
    return sigm_s,pix_s

def weighterror1(St,Wt,multiT=True,nprocf=6):
    nly,nlx,ns=Wt.shape
    out=np.zeros([nly,nlx,ns])
    if multiT:
        nproc=nprocf
        with ThreadPool(nproc) as pool:
            args=[]
            for npros in range(0, nproc):
                val=int(nlx/nproc)
                a1=val*npros
                if npros < nproc-1:
                    a2=val*(npros+1)
                else:
                    a2=nlx
                Wgt=Wt[:,a1:a2,:]     
                args.extend([(St,Wgt,nly,ns,a1,a2)])                    
            result_l = pool.map(kernel.task_wrappercov1, args)
    else:
        nproc=1
        npros=0
        result_l=[]
        args=(St,Wt,nly,ns,0,nlx)
        result_l.extend([kernel.task_wrappercov1(args)])
    for npros in range(0, nproc):
        result=result_l[npros]
        val=int(nlx/nproc)
        a1=val*npros
        if npros < nproc-1:
            a2=val*(npros+1)
        else:
            a2=nlx
        out[:,a1:a2,:]=result[0]
    return out

def weighterror2(St,Wt,multiT=True,nprocf=6,verbose=True,matf=True):
    #This function will propagate all the error and generate the full covariance matrix 
    #St is the error matrix frim all the fibres
    #Wt are the weigths for the 2d image interpolation
    #if matf eq False the output shape will be [nx,ny,nx,ny], if True the shape will be [nx*ny,nx*ny] 
    out=weighterror1(St,Wt,multiT=multiT,nprocf=nprocf)
    nly,nlx,ns=Wt.shape
    if matf:
        outf=np.zeros([nly*nlx,nly*nlx])
        dist=np.zeros([nly*nlx,nly*nlx])
    else:
        outf=np.zeros([nly,nlx,nly,nlx])
        dist=np.zeros([nly,nlx,nly,nlx]) 
    indexT=np.array([(k,0) for k in range(nly)])
    Dq=pdist(indexT, metric='euclidean')
    if verbose:
        pbar=tqdm(total=nlx)
    for i in range(0, nlx):
        Wgt=Wt[:,i,:]
        if multiT:
            nproc=nprocf
            with ThreadPool(nproc) as pool:
                args=[]
                for npros in range(0, nproc):
                    val=int(nlx/nproc)
                    a1=val*npros
                    if npros < nproc-1:
                        a2=val*(npros+1)
                    else:
                        a2=nlx
                    outt=out[:,a1:a2,:] 
                    args.extend([(St,Wgt,outt,Dq,nly,i,a1,a2)])                   
                result_l = pool.map(kernel.task_wrappercov2, args)
        else:
            nproc=1
            npros=0
            result_l=[]
            args=(St,Wgt,out,Dq,nly,i,0,nlx)
            result_l.extend([kernel.task_wrappercov2(args)])
        for npros in range(0, nproc):
            result=result_l[npros]
            val=int(nlx/nproc)
            a1=val*npros
            if npros < nproc-1:
                a2=val*(npros+1)
            else:
                a2=nlx
            if matf:
                b1=i*nly
                b2=(i+1)*nly
                for j in range(a1,a2):
                    c1=j*nly
                    c2=(j+1)*nly
                    outf[b1:b2,c1:c2]=result[0][:,:,j-a1]
                    dist[b1:b2,c1:c2]=result[1][:,:,j-a1]
            else:
                outf[:,i,:,a1:a2]=result[0]
                dist[:,i,:,a1:a2]=result[1]
        if verbose:        
            pbar.update(1)
    if verbose:
        pbar.close()    
    return outf,dist

def correlation_matrix(out):
    nx,ny=out.shape
    errf0=np.ones([nx])
    outf0=np.zeros([nx,nx])
    for i in range(0, nx):
        if out[i,i] > 0:
                errf0[i]=out[i,i]
    for i in range(0, nx):
        for j in range(0, nx):
            outf0[i,j]=out[i,j]/np.sqrt(errf0[i]*errf0[j])
    return outf0   

def get_error(out,Wg1):
    nx,ny,ns=Wg1.shape
    errt=np.zeros([nx,ny])
    for i in range(0, ny):
        b1=i*nx
        b2=(i+1)*nx
        for k in range(b1,b2):
            errt[k-b1,i]=out1[(k+0) % (nx*ny),(k+0) % (nx*ny)]
    return errt
