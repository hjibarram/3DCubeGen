import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
import CubeGen.tools.tools as tools
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy import units as u
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
import os.path as ptt
import json
import copy
import shutil

def extintion_c(wave,dir_tem='data',basename='extintion_curve'):
    f=open(dir_tem+'/'+basename+'.txt','r')
    x1=[]
    y1=[]
    for line in f:
        data=line.replace('\n','').split(' ')
        data=list(filter(None,data))
        x1.extend([float(data[0])])
        y1.extend([float(data[1])])
    y1=np.array(y1)
    x1=np.array(x1)
    Kl=interp1d(x1,y1,bounds_error=False,fill_value=0.)(wave)
    return Kl

def id_str(id,n_z=2):
    id=int(float(id))
    if n_z < 2 or n_z > 8:
        n_z=2
    if n_z == 2:
        if id < 10:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 3:
        if id < 10:
            idt='00'+str(id)
        elif id < 100:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 4:
        if id < 10:
            idt='000'+str(id)
        elif id < 100:
            idt='00'+str(id)
        elif id < 1000:
            idt='0'+str(id)    
        else:
            idt=str(id)        
    return idt

def megarafiber_pos(hdr,verbose=False,astmet=True):
    nfib=hdr['NFIBERS']
    psc=hdr['PSCALE']
    x_pos=np.zeros(nfib)
    y_pos=np.zeros(nfib)
    fib_a=np.zeros(nfib)
    fib_b=np.zeros(nfib)
    fib_id=np.zeros(nfib)
    for i in range(0, nfib):
        x_pos[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_X']
        y_pos[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_Y']
        fib_a[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_D']
        fib_b[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_B']
        fib_id[i]=i+1
    nt=np.where((np.abs(x_pos) < 10) & (np.abs(y_pos) < 10))
    x_posf=x_pos[nt]
    y_posf=y_pos[nt]
    fib_idt=fib_id[nt]
    nt=np.where((np.abs(x_pos) > 10) & (np.abs(y_pos) > 10))
    fib_ids=fib_id[nt]
    if astmet:
        try:
            pc11=hdr['PC1_1']
            pc12=hdr['PC1_2']
            pc21=hdr['PC2_1']
            pc22=hdr['PC2_2']
            x_ifu=(x_posf*pc11+y_posf*pc21)*psc+hdr['CRVAL1']*3600.0    
            y_ifu=(x_posf*pc12+y_posf*pc22)*psc+hdr['CRVAL2']*3600.0
        except:
            x_ifu=x_posf*psc+hdr['CRVAL1']*3600.0
            y_ifu=y_posf*psc+hdr['CRVAL2']*3600.0 
    else:
        x_ifu=x_posf*psc
        y_ifu=y_posf*psc
    if verbose:
        import matplotlib.pyplot as plt
        plt.plot(x_ifu/3600.,y_ifu/3600.,'o')
        plt.show() 
    return x_ifu,y_ifu,fib_idt,fib_ids

def read_standar(path_data='data',stdar_t='Feige32',stdT='',fergs=True):
    wav_i=[]
    res_i=[]
    file=path_data+'/'+stdar_t+stdT+'.dat'
    f=open(file,'r')
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=list(filter(None,data))
            wav_i.extend([float(data[0])])
            if fergs:
                res_i.extend([float(data[1])*1e-16])
            else:
                res_i.extend([float(data[1])])
    f.close()
    wav_i=np.array(wav_i)
    res_i=np.array(res_i)
    return wav_i,res_i

def extract_cube(wave0,s,r,xo=0,yo=0,stdar_t='Feige32',path_ifu='',vph='B',dpix=0.35):
    file=path_ifu+'/SPSTD_'+stdar_t+'_'+vph+'.fits.gz'
    [flux, hdrT]=fits.getdata(file, 0, header=True)
    nz,nx,ny=flux.shape
    rad=np.zeros([nx,ny])
    for i in range(0, nx):
        for j in range(0,ny):
            rt=np.sqrt((i-xo)**2.0+(j-yo)**2.0)*dpix
            if rt <= r:
                rad[i,j]=1.0
    fluxf=np.zeros(nz)
    nt=np.where(rad == 1.0)
    for z in range(0, nz):
        fluxf0=flux[z,:,:]
        fluxf[z]=np.nansum(fluxf0[nt])
    flux=fluxf
    l0=hdrT['CRVAL3']
    dl=hdrT['CDELT3']
    wave=l0+np.arange(len(flux))*dl
    s0=np.copy(s)
    s=interp1d(wave0,s0,bounds_error=False,fill_value=0.)(wave)
    flux=flux*1e-16#/1.501045#/0.172
    return flux,wave,s

def get_rssflux_sens(r,xo=0,yo=0,path_block9='',path_sensfits='',path_data='data',vph='B',phase=0,dpix=0.35):
    file=path_sensfits+'master_sensitivity.fits'
    if phase >= 1:
        file=path_data+'/master_sensitivity'+vph+'_n.fits'
    [s, hdr2]=fits.getdata(file, 0, header=True)
    file=path_block9+'final_rss.fits'
    [flux, hdr0]=fits.getdata(file, 0, header=True)
    l0=hdr0['CRVAL1']
    dl=hdr0['CDELT1']
    Xa=hdr0['AIRMASS']
    [flux0, hdr1]=fits.getdata(file, 1, header=True)
    x_ifu,y_ifu,fib_idt,fib_ids=megarafiber_pos(hdr1,astmet=False)
    xs=xo*dpix-12.5/2#7 #-npix/2*dpix IFU
    ys=yo*dpix-11.3/2#7 #-npix/2*dpix IFU
    nt=np.where(np.sqrt((x_ifu-xs)**2.0+(y_ifu-ys)**2.0) <= r)
    fib_idt=fib_idt[nt]
    fluxf=0
    for i in range(0, len(fib_idt)):
        fluxf=fluxf+flux[int(fib_idt[i])-1,:]
    flux=fluxf
    wave=l0+np.arange(len(flux))*dl
    return wave,flux,Xa,s,hdr2


def high_res_std(wav_i,res_i,wave,flux,flux1,path_data='',stdar_t='Feige32',vph='B'):
	#Generate high resolution standard star spectra
    #First line
    la1=4800#7550#6530#6884
    la2=4900#7619#6590#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    flutt1=res_i[nt1]
    flutt2=res_i[nt2]
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    flutt0=flux1[nt0]
    #flutt0[len(flutt0)-2:len(flutt0)]=4.23e-11#np.mean(flutt0[0:3])#esto es valido solo para este objeto
    wavt0=wave[nt0]
    res_i=np.append(np.append(flutt1,flutt0),flutt2)
    wav_i=np.append(np.append(wavt1,wavt0),wavt2)
    #Second line
    la1=6530
    la2=6590
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    flutt1=res_i[nt1]
    flutt2=res_i[nt2]
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    flutt0=flux1[nt0]
    wavt0=wave[nt0]
    res_i=np.append(np.append(flutt1,flutt0),flutt2)
    wav_i=np.append(np.append(wavt1,wavt0),wavt2)
    #Third line
    la1=7470#6884
    la2=7500#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    flutt1=res_i[nt1]
    flutt2=res_i[nt2]
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    flutt0=flux1[nt0]
    wavt0=wave[nt0]
    res_i=np.append(np.append(flutt1,flutt0),flutt2)
    wav_i=np.append(np.append(wavt1,wavt0),wavt2)
    #First Atmospheric bands 
    la1=6860#7619#6860#6884
    la2=6950#7720#6950#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    wavt0=wave[nt0]
    wav_it=np.append(np.append(wavt1,wavt0),wavt2)
    nt3=np.where((wav_i <= la1) | (wav_i >= la2))
    flutt3=res_i[nt3]
    wavt3=wav_i[nt3]
    fluxft=interp1d(wavt3,flutt3,bounds_error=False,fill_value=0.)(wav_it)
    res_i=fluxft
    wav_i=wav_it
    #Second Atmospheric bands 
    la1=6950#8100#6950#7150#6884
    la2=7380#8350#7380#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    wavt0=wave[nt0]
    wav_it=np.append(np.append(wavt1,wavt0),wavt2)
    nt3=np.where((wav_i <= la1) | (wav_i >= la2))
    flutt3=res_i[nt3]
    wavt3=wav_i[nt3]
    fluxft=interp1d(wavt3,flutt3,bounds_error=False,fill_value=0.)(wav_it)
    res_i=fluxft
    wav_i=wav_it
    #Third Atmospheric bands Add commentMore actions
    la1=7550#8100#6950#7150#6884
    la2=7800#8350#7380#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    wavt0=wave[nt0]
    wav_it=np.append(np.append(wavt1,wavt0),wavt2)
    nt3=np.where((wav_i <= la1) | (wav_i >= la2))
    flutt3=res_i[nt3]
    wavt3=wav_i[nt3]
    fluxft=interp1d(wavt3,flutt3,bounds_error=False,fill_value=0.)(wav_it)
    res_i=fluxft
    wav_i=wav_it
    #Forth Atmospheric bands 
    la1=8100#6950#7150#6884
    la2=8350#7380#7120
    nt1=np.where((wav_i <= la1))
    nt2=np.where((wav_i >= la2))
    wavt1=wav_i[nt1]
    wavt2=wav_i[nt2]
    nt0=np.where((wave > la1) & (wave < la2))
    wavt0=wave[nt0]
    wav_it=np.append(np.append(wavt1,wavt0),wavt2)
    nt3=np.where((wav_i <= la1) | (wav_i >= la2))
    flutt3=res_i[nt3]
    wavt3=wav_i[nt3]
    fluxft=interp1d(wavt3,flutt3,bounds_error=False,fill_value=0.)(wav_it)
    res_i=fluxft
    wav_i=wav_it
    f=open(path_data+'/'+stdar_t+'_hr'+vph+'.dat','w')
    for i in range(0, len(wav_i)):
        f.write(str(wav_i[i])+' '+str(res_i[i])+' \n')
    f.close()
    return wav_i,res_i
    #end high resolution standard star spectra


def gen_sensf(wav_i,wave,flux,res_f,at_ext,minwave,maxwave,hdr2,path_data='data',vph='B',scalefact=1.0):
    la1=6550#7594#6550
    la2=6590#7619#6590
    lb1=6860#8100#6860
    lb2=6910#8300#6910#6890
    nt=np.where((wav_i >= np.min(wave)) & (wav_i <= np.max(wave)))
    fluxt=interp1d(wave,flux,bounds_error=False,fill_value=0.)(wav_i[nt])
    fluxt1=interp1d(wav_i[nt],fluxt,bounds_error=False,fill_value=0.)(wave)
    res_f=res_f*at_ext
    resp=fluxt1/res_f
    nt=np.where((wave >= minwave-20) & (wave <= maxwave+20))
    resp_t=resp[nt]
    resp_t0=resp[nt]
    wave_t=wave[nt]

    nt0=np.where((wave_t >= la1) & (wave_t <= la2))
    resp_t[nt0]=0
    nt1=np.where((wave_t >= lb1) & (wave_t <= lb2))
    resp_t[nt1]=0
    nt2=np.where(resp_t > 0)
    resp_tt=resp_t[nt2]
    wave_tt=wave_t[nt2]
    resp_t=interp1d(wave_tt,resp_tt,bounds_error=False,fill_value=0.)(wave_t)
    resp_t1=tools.median_a(resp_t,lw=60)
    resp_t1[nt0]=0
    resp_t1[nt1]=0
    resp_t1=resp_t1[nt2]
    resp_t1=interp1d(wave_tt,resp_t1,bounds_error=False,fill_value=0.)(wave_t)

    nt3=np.where((wave_t >= 8100) & (wave_t <= 8400))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    nt3=np.where((wave_t >= 6860) & (wave_t <= 7400))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    #nt3=np.where((wave_t >= 7617) & (wave_t <= 7720))#6884-7400
    nt3=np.where((wave_t >= 6860) & (wave_t <= 7400))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    nt3=np.where((wave_t >= 7500) & (wave_t <= 7780))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    resp[nt]=resp_t1#model
    s=resp
    fits.writeto(path_data+'/master_sensitivity'+vph+'_n.fits', resp*scalefact, hdr2,overwrite=True)#*0.2*0.555#*0.3336#*0.864
    flux=flux/s
    return flux,s,res_f

def meg_spectra_outs(wave,flux,res_f,minwave,maxwave,vph='B',pathwork=''):
    nt=np.where((wave > minwave) & np.isfinite(flux))# & (s > 0)) 4332
    max_val1=np.nanmax(flux[nt])*1.1
    factor=np.nanmean(flux[nt]/res_f[nt])
    fig, ax = plt.subplots(figsize=(6.5,5.5))
    ax.set_ylim(0,max_val1)#/100.0)
    ax.set_xlim(minwave,maxwave)#4332,5230)#7250,8700)#6109,7399)#3500,10100)  7550,7720)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel(r'Flux $[erg/s/cm^2/\AA]$',fontsize=14)
    plt.plot(wave,flux,alpha=0.4,color='grey')
    #plt.plot(wave,s,alpha=1.0,color='black')#spectro final
    plt.plot(wave,res_f,color='blue',alpha=0.8)#spectro original
    #plt.plot(wave_t,resp_t1,alpha=0.9,color='red')
    #plt.plot(wave_t,resp_t,alpha=0.9,color='green')
    fig.tight_layout()
    plt.savefig(pathwork+'spec'+vph+'.jpg')#,dpi=1000)
    plt.close()
    return factor

def calib_spec(phase=0,scfact=1.0,xo=0,yo=0,vph='B',idp='9',stdar_t='Feige32',verbose=True,path_data='data',pathwork='',path_ifu='ifu',r=4.0,dpix=0.35):
    #r is Aperture radius
    path_block9=pathwork+'obsid'+idp+vph+'_results/'
    path_sensfits=pathwork+'obsid'+idp+vph+'_results/'
    if phase == 0:
        cube=False
        gen_hr=False
        gensensf=True
        scalefact=1.0
        fergs=True
        stdT=''
    if phase == 1:
        cube=False
        gen_hr=True
        gensensf=False
        scalefact=1.0
        fergs=True
        stdT=''
    if phase == 2:
        cube=False
        gen_hr=False
        gensensf=True
        scalefact=1.0
        fergs=False
        stdT='_hr'+vph
    if phase == 3 or phase == 5:
        cube=True
        gen_hr=False
        gensensf=False
        scalefact=1.0
        fergs=False
        stdT='_hr'+vph
    if phase == 4:
        cube=False
        gen_hr=False
        gensensf=True
        scalefact=scfact
        fergs=False
        stdT='_hr'+vph
    wav_i,res_i=read_standar(path_data=path_data,stdar_t=stdar_t,stdT=stdT,fergs=fergs)
    wave,flux,Xa,s,hdr2=get_rssflux_sens(r,xo=xo,yo=yo,path_block9=path_block9,path_sensfits=path_sensfits,path_data=path_data,vph=vph,phase=phase,dpix=dpix)
    if cube == True:
        flux,wave,s=extract_cube(wave,s,r,xo=xo,yo=yo,stdar_t=stdar_t,path_ifu=path_ifu,vph=vph,dpix=dpix)
    maxwave=np.round(np.nanmax(wave)-70)
    minwave=np.round(np.nanmin(wave)+70)
    Kvl=extintion_c(wave,dir_tem=path_data)
    at_ext=10.0**(-0.4*Xa*Kvl)
    flux1=flux/s/at_ext
    if gen_hr == True:
        wav_i,res_i=high_res_std(wav_i,res_i,wave,flux,flux1,path_data=path_data,stdar_t=stdar_t,vph=vph)
    res_f=interp1d(wav_i,res_i,bounds_error=False,fill_value=0.)(wave)
    if gensensf == True:
        flux,s,res_f=gen_sensf(wav_i,wave,flux,res_f,at_ext,minwave,maxwave,hdr2,path_data=path_data,vph=vph,scalefact=scalefact)
    factor=meg_spectra_outs(wave,flux,res_f,minwave,maxwave,vph=vph,pathwork=pathwork)
    if verbose:
        print('The scale factor is:',factor)
    return factor


def lst_c(mjd=56000.0,tai=np.inf,lng=-105.820417):
    if np.isfinite(tai) == True:
        jd=2400000.5+tai/(24.0*3600.0)
    else:
        jd=mjd+2400000.5
    c = [280.46061837,360.98564736629,0.000387933,38710000.0]
    jd2000=2451545.0
    t0 = jd - jd2000
    t = t0/36525.0
    theta=c[0]+(c[1]*t0)+t**2.0*(c[2]-t/c[3])
    
    lst=(theta+lng+360.0)/15.0
    if lst < 0.0:
        lst = 24.0+(lst % 24.0)
    lst = lst % 24.0
    return lst

def paralactic_angle(ha,dec=0.0,phi=32.78):
    p=np.arctan2(np.sin(ha/24.0*360.0*np.pi/180.0),(np.tan(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)-np.sin(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)))*180.0/np.pi
    p=p % 360.0
    return p    

def airmas(ha,dec=0.0,phi=32.78):
    cosz=np.sin(phi*np.pi/180.0)*np.sin(dec*np.pi/180.0)+np.cos(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)
    M=(1.002432*cosz**2.0+0.148386*cosz+0.0096467)/(cosz**3.0+0.149864*cosz**2.0+0.0102963*cosz+0.000303978)# Airmas formula from Young 1994
    return M
   
def refrac_dif(wave,ha,dec=0.0,phi=32.78,lo=5070.0,T=7.0,P=600.0,f=8.0,vapor=False):
    cosz=np.sin(phi*np.pi/180.0)*np.sin(dec*np.pi/180.0)+np.cos(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)
    R=64.328+29498.1/(146.0-(1.0/(wave*1e-4))**2.0)+255.4/(41.0-(1.0/(wave*1e-4))**2.0)
    R=R*P*(1.0+P*(1.049-0.0157*T)*1e-6)/(720.883*(1.0+0.003661*T))
    R1=64.328+29498.1/(146.0-(1.0/(lo*1e-4))**2.0)+255.4/(41.0-(1.0/(lo*1e-4))**2.0)
    R1=R1*P*(1.0+P*(1.049-0.0157*T)*1e-6)/(720.883*(1.0+0.003661*T))
    if vapor == True:
        R=R*(0.0624-0.000680/(wave*1e-4)**2.0)/(1.0+0.003661*T)*f
        R1=R1*(0.0624-0.000680/(lo*1e-4)**2.0)/(1.0+0.003661*T)*f
    R=R-R1
    z=np.arccos(cosz)
    R=R*np.tan(z)
    R=R/1e6*3600*180/np.pi
    return R


def wavelength_virus(header):
    ref_pixel = header['CRPIX1']
    coord_ref_pixel = header['CRVAL1']
    wave_pixel = header['CDELT1']
    npix=header['NAXIS1']
    wstart = coord_ref_pixel - ( (ref_pixel-1)*wave_pixel)
    wave = np.array([wstart + i*wave_pixel for i in range(npix)])
    return wave,ref_pixel,coord_ref_pixel,wave_pixel

def get_adr(hdr,wave,lo=5000.0,repss=True):
    mjd=hdr["MJD-OBS"]
    if repss:
        lng=hdr["LONGITUD"].replace("+","-")
    else:
        lng=hdr["LONGITUD"]
    lat=hdr["LATITUDE"]
    hig=hdr["HEIGHT"]
    Presure=hdr["PRESSURE"]
    Temp=hdr["TAMBIENT"]
    Tai=hdr["DATE-OBS"]
    pos_obs = EarthLocation(lng,lat,hig*u.m)
    long=pos_obs.lon.deg
    lati=pos_obs.lat.deg
    min_wave=np.nanmin(wave)
    max_wave=np.nanmax(wave)    
    #time=Time(Tai,format="fits",scale="tai")
    #tai=time.to_value('unix_tai', 'long')
    #tai=time.to_value('mjd', 'long')*(24.0*3600.0
    dec_0=hdr["DECDEG"]
    ra_0=hdr["RADEG"]
    ha=lst_c(lng=long,mjd=mjd)-ra_0/180.0*12.0
    R=refrac_dif(wave,ha,dec=dec_0,phi=lati,lo=lo,T=Temp,P=Presure,f=8.0)
    pa=paralactic_angle(ha,dec=dec_0,phi=lati)
    R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
    return R2,R

def get_radvel(hdr,repss=True):
    import astropy.units as u
    mjd=hdr["MJD-OBS"]
    if repss:
        lng=hdr["LONGITUD"].replace("+","-")
    else:
        lng=hdr["LONGITUD"]
    lat=hdr["LATITUDE"]
    hig=hdr["HEIGHT"]
    Tai=hdr["DATE-OBS"]
    dec_0=hdr["DECDEG"]
    ra_0=hdr["RADEG"]
    sky = SkyCoord(ra_0,dec_0,unit="deg",frame="fk5")
    pos_obs = EarthLocation(lng,lat,hig*u.m)
    time = Time(mjd, format='mjd', scale='utc')
    rv=sky.radial_velocity_correction(obstime=time, location=pos_obs).value 
    rv=rv/299792458.0
    return rv


def evaluate_2dPSF(pf_map,model=True,sig=2,plotview=False,centroid=False):
    nx,ny=pf_map.shape
    if sig == 0:
        pf_map_c=pf_map
    else:
        PSF=Gaussian2DKernel(x_stddev=sig,y_stddev=sig)
        pf_map_c=convolve(pf_map, PSF)
    min_in=np.unravel_index(np.nanargmax(pf_map_c), (nx,ny))
    At=np.nanmax(pf_map)
    x_t=np.arange(ny)-min_in[1]
    y_t=np.arange(nx)-min_in[0]
    x_t=np.array([x_t]*nx)
    y_t=np.array([y_t]*ny).T
    n_ds=100
    n_dx=50
    n_dy=50
    ds_t=np.arange(n_ds)/float(n_ds)*(4.0-1.0)+1.0
    dx=np.arange(n_dx)/float(n_dx)*10-5.0
    dy=np.arange(n_dy)/float(n_dy)*10-5.0
    chi2_m=np.zeros([n_ds,n_dx,n_dy])
    for j in range(0, n_ds):
        for k in range(0, n_dx):
            for w in range(0, n_dy):
                spec_t=np.exp(-0.5*((((x_t-dx[k])/ds_t[j])**2.0)+((y_t-dy[w])/ds_t[j])**2.0))*At
                chi2_m[j,k,w]=np.nansum((pf_map-spec_t)**2.0)
    min_int=np.unravel_index(np.nanargmin(chi2_m), (n_ds,n_dx,n_dy))
    ds_m=ds_t[min_int[0]]
    dx_m=dx[min_int[1]]
    dy_m=dy[min_int[2]]
    if plotview:
        cm=plt.cm.get_cmap('jet')
        lev=np.sqrt(np.arange(0.0,10.0,1.5)+0.008)/np.sqrt(10.008)
        fig, ax = plt.subplots(figsize=(6.8*1.1,5.5*1.2))
        #ict=plt.imshow(np.log10(pf_map/At),cmap=cm) 
        ict=plt.imshow((pf_map_c/At),cmap=cm) 
        cbar=plt.colorbar(ict)
        ics=plt.contour(pf_map/At,lev,colors='k',linewidths=1)            
        cbar.set_label(r"Relative Density")
        fig.tight_layout()
        plt.show()
    if model:
        spec_t=np.exp(-0.5*((((x_t-dx_m)/ds_m)**2.0)+((y_t-dy_m)/ds_m)**2.0))*At
        if plotview:
            fig, ax = plt.subplots(figsize=(6.8*1.1,5.5*1.2))
            ict=plt.imshow(spec_t/At,cmap=cm) 
            cbar=plt.colorbar(ict)
            ics=plt.contour(spec_t/At,lev,colors='k',linewidths=1)
            ics=plt.contour(pf_map/At,lev,colors='red',linewidths=1)            
            cbar.set_label(r"Relative Density")
            fig.tight_layout()
            plt.show()
            fig, ax = plt.subplots(figsize=(6.8*1.1,5.5*1.2))
            ict=plt.imshow((pf_map-spec_t)/At,cmap=cm) 
            cbar=plt.colorbar(ict)
            ics=plt.contour((pf_map-spec_t)/At,lev,colors='k',linewidths=1)
            cbar.set_label(r"Relative Density")
            fig.tight_layout()
            plt.show()
        dx_m=dx_m+min_in[1]#0
        dy_m=dy_m+min_in[0]#1
        psf=ds_m*2.0*np.sqrt(2.0*np.log10(2.0))
        if centroid:
            xo=dx_m
            yo=dy_m
            return xo,yo
        else:
            return dx_m,dy_m,ds_m,psf,spec_t
    else:    
        dx_m=dx_m+min_in[1]
        dy_m=dy_m+min_in[0]
        psf=ds_m*2.0*np.sqrt(2.0*np.log10(2.0))
        if centroid:
            xo=dx_m
            yo=dy_m
            return xo,yo
        else:
            return dx_m,dy_m,ds_m,psf

def read_obs(name,path=''):
    """
    Reads an observation file and returns the data as a numpy array.
    Parameters:
    name (str): The name of the file to read.
    path (str): The directory path where the file is located.

    Returns:
    np.ndarray: The data read from the file.
    """
    file=path+name
    f=open(file,'r')
    dic={}
    ct=0
    img=1000
    for line in f:
        ct+=1
        if 'IMAGE' in line:
            data=line.replace('\n','').split(' ')
            data=list(filter(None, data)) # Remove empty strings
            data.extend(['Type'])
            nh=len(data)
            for it in range(0, nh):
                dic.update({data[it]:[]})
            head=data
            img=ct
            typ='IMAGE'
        if 'Flat' in line and img != 1000:
            typ='Flat'
            img=ct
        if 'Bias' in line and img != 1000 and not '.fits' in line:
            typ='Bias'
            img=ct
        if 'Arcs' in line and img != 1000:
            typ='Arcs'
            img=ct
        if (('Spectrophotometric Standard' in line) or ('Spectrophotometric standard' in line) or ('Stds' in line)) and img != 1000:    
            typ='SPTD'
            img=ct
        if ct >=  img+2:
            data=line.replace('\n','').split(' ')
            data=list(filter(None,data))
            if len(data)+1 == nh:
                for it in range(0, nh-1):
                    try:
                        if it == 3:
                            val=data[it].replace(' ','')
                        else:
                            val=float(data[it])
                    except:
                        val=data[it].replace(' ','')
                    dic[head[it]].extend([val])
                dic[head[it+1]].extend([typ])
    f.close()
    return dic,head
            
def create_obsBIAS(data,redux_path='',ob_path=''):
    """
    Creates a bias observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where(np.array(data['Type']) == 'Bias')[0]
    if len(nt) == 0:
        print('No Bias observations found')
        return
    dat=np.array(data['IMAGE'])
    filelists=dat[nt]
    file=redux_path+'/obsresult-1.yaml'
    f=open(file,'w')
    f.write('id: 1\n')
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraBiasImage\n')
    f.write('images:\n')
    for i in range(0, len(nt)):
        line='  - '+filelists[i]+'\n'
        f.write(line)
        if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
            if ptt.exists(ob_path+'bias/'+filelists[i]+'.gz') == True:
                call='gunzip -c '+ob_path+'bias/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
            else:
                call='cp '+ob_path+'bias/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
            tools.sycall(call)
    f.close()
    
def get_vph(data):
    """
    Returns the VPH (Variable Passband Holography) value from the observation file data.
    data (dic): Dictionary from the observation file data

    Returns:
    str: The VPH value.
    """
    nt=np.where(np.array(data['Type']) == 'IMAGE')[0]
    dat=np.array(data['VPH'])
    vph=np.unique(dat[nt])
    return vph

def get_std(data):
    """
    Returns the standard star name from the observation file data.
    Parameters:
    data (dic): Dictionary from the observation file data.

    Returns:
    str: The standard star name.
    """
    nt=np.where(np.array(data['Type']) == 'SPTD')[0]
    dat=np.array(data['OBJECT'])
    std=dat[nt]
    if len(std) > 0:
        std=np.unique(std)
        for i in range(0, len(std)):
            std[i]=std[i].replace('SPSTD_','')
        return std
    else:
        return None

def get_exptime(data):
    """
    Returns the exposure time from the observation file data.
    Parameters:
    data (dic): Dictionary from the observation file data.

    Returns:
    float: The exposure time in seconds.
    """
    nt=np.where(np.array(data['Type']) == 'IMAGE')[0]
    dat=np.array(data['EXPTIME'])
    dat1=np.array(data['VPH'])
    exptime=dat[nt]
    vph=dat1[nt]
    if len(exptime) > 0:
        vphu=np.unique(vph)
        exptimeT= np.zeros(len(vphu))
        for j in range(0, len(vphu)):
            for i in range(0, len(exptime)):
                if vph[i] == vphu[j]:
                    exptimeT[j]=exptime[i]+exptimeT[j]
        return vphu,exptimeT
    else:
        return None


def get_obj(data):
    """
    Returns the Object name from the observation file data.
    Parameters:
    data (dic): Dictionary from the observation file data.

    Returns:
    str: The Object name.
    """
    nt=np.where(np.array(data['Type']) == 'IMAGE')[0]
    dat=np.array(data['OBJECT'])
    obj=dat[nt]
    if len(obj) > 0:
        obj=np.unique(obj)
        for i in range(0, len(obj)):
            obj[i]=obj[i].replace('SPSTD_','')
        return obj
    else:
        return None

def create_obsTraceMap(data,vph,redux_path='',ob_path=''):
    """
    Creates a TraceMap observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'Flat'))[0]
    if len(nt) == 0:
        print('No Flat observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    file=redux_path+'/obsresult-2VPH.yaml'.replace('VPH',vph)
    f=open(file,'w')
    f.write('id: 2VPH\n'.replace('VPH',vph))
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraTraceMap\n')
    f.write('images:\n')
    for i in range(0, len(filelists)):
        if vphs[i] == vph:
            line='  - '+filelists[i]+'\n'
            f.write(line)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                if ptt.exists(ob_path+'flat/'+filelists[i]+'.gz') == True:
                    call='gunzip -c '+ob_path+'flat/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                else:
                    call='cp '+ob_path+'flat/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                tools.sycall(call)
    f.close()

def create_obsFiberFlatImage(data,vph,redux_path='',ob_path=''):
    """
    Creates a FiberFlatImage observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'Flat'))[0]
    if len(nt) == 0:
        print('No Flat observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    file=redux_path+'/obsresult-6VPH.yaml'.replace('VPH',vph)
    f=open(file,'w')
    f.write('id: 6VPH\n'.replace('VPH',vph))
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraFiberFlatImage\n')
    f.write('images:\n')
    for i in range(0, len(filelists)):
        if vphs[i] == vph:
            line='  - '+filelists[i]+'\n'
            f.write(line)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                if ptt.exists(ob_path+'flat/'+filelists[i]+'.gz') == True:
                    call='gunzip -c '+ob_path+'flat/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                else:
                    call='cp '+ob_path+'flat/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                tools.sycall(call)
    f.close()

def create_obsModelMap(data,vph,redux_path='',ob_path=''):
    """
    Creates a ModelMap observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'Flat'))[0]
    if len(nt) == 0:
        print('No Flat observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    file=redux_path+'/obsresult-8VPH.yaml'.replace('VPH',vph)
    f=open(file,'w')
    f.write('id: 8VPH\n'.replace('VPH',vph))
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraModelMap\n')
    f.write('images:\n')
    for i in range(0, len(filelists)):
        if vphs[i] == vph:
            line='  - '+filelists[i]+'\n'
            f.write(line)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                if ptt.exists(ob_path+'flat/'+filelists[i]+'.gz') == True:
                    call='gunzip -c '+ob_path+'flat/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                else:
                    call='cp '+ob_path+'flat/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                tools.sycall(call)
    f.close()

def create_obsArcCalibration(data,vph,redux_path='',ob_path='',forceAr=False):
    """
    Creates a ArcCalibration observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'Arcs'))[0]
    if len(nt) == 0:
        print('No Arcs observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    dat2=np.array(data['OSFILT'])
    dat3=np.array(data['ThNe1on'])
    dat4=np.array(data['ThAr1on'])
    dat5=np.array(data['ThAr3on'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filt=dat2[nt]
    ThNe=dat3[nt]
    ThAr=dat4[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    filt=filt[uin]
    ThNe=ThNe[uin]
    ThAr=ThAr[uin]
    ThAr3=dat5[nt]
    file=redux_path+'/obsresult-3VPH.yaml'.replace('VPH',vph)
    f=open(file,'w')
    f.write('id: 3VPH\n'.replace('VPH',vph))
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraArcCalibration\n')
    f.write('images:\n')
    for i in range(0, len(filelists)):
        if filt[i] == 'BLUE' and ThAr[i] == 1:
            svt=True
        elif filt[i] == 'RED' and ThNe[i] == 1:
            svt=True
        else:
            svt=False
        if forceAr:
            if filt[i] == 'BLUE' and ThAr[i] == 1:
                svt=True
            elif filt[i] == 'RED' and ThAr3[i] == 1:
                svt=True
            else:
                svt=False    
        if vphs[i] == vph and svt:
            line='  - '+filelists[i]+'\n'
            f.write(line)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                if ptt.exists(ob_path+'arc/'+filelists[i]+'.gz') == True:
                    call='gunzip -c '+ob_path+'arc/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                else:
                    call='cp '+ob_path+'arc/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                tools.sycall(call)
    f.close()


def create_obsLcbStdStar(data,vph,redux_path='',ob_path=''):
    """
    Creates a LcbStdStar observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'SPTD'))[0]
    if len(nt) == 0:
        print('No SPTD observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    dat2=np.array(data['OBJECT'])
    filelists=dat[nt]
    vphs=dat1[nt]
    stdt=dat2[nt]
    ustdt=np.unique(stdt)
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    stdt=stdt[uin]
    idst=['9','10','14','15']
    for j in range(0, len(ustdt)):
        file=redux_path+'/obsresult-NVPH.yaml'.replace('VPH',vph).replace('N',idst[j])
        f=open(file,'w')
        f.write('id: NVPH\n'.replace('VPH',vph).replace('N',idst[j]))
        f.write('instrument: MEGARA\n')
        f.write('mode: MegaraLcbStdStar\n')
        f.write('images:\n')
        for i in range(0, len(filelists)):
            if stdt[i] == ustdt[j]:
                svt=True
            else:
                svt=False
            if vphs[i] == vph and svt:
                line='  - '+filelists[i]+'\n'
                f.write(line)
                if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                    if ptt.exists(ob_path+'stds/'+filelists[i]+'.gz') == True:
                        call='gunzip -c '+ob_path+'stds/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                    else:
                        call='cp '+ob_path+'stds/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                    tools.sycall(call)
        f.close()

def create_obsLcbImageStd(data,vph,redux_path='',ob_path=''):
    """
    Creates a LcbImageStd observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'SPTD'))[0]
    if len(nt) == 0:
        print('No SPTD observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    dat2=np.array(data['OBJECT'])
    filelists=dat[nt]
    vphs=dat1[nt]
    stdt=dat2[nt]
    ustdt=np.unique(stdt)
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    stdt=stdt[uin]
    idst=['11','13','16','17']
    for j in range(0, len(ustdt)):
        file=redux_path+'/obsresult-NVPH.yaml'.replace('VPH',vph).replace('N',idst[j])
        f=open(file,'w')
        f.write('id: NVPH\n'.replace('VPH',vph).replace('N',idst[j]))
        f.write('instrument: MEGARA\n')
        f.write('mode: MegaraLcbImage\n')
        f.write('images:\n')
        for i in range(0, len(filelists)):
            if stdt[i] == ustdt[j]:
                svt=True
            else:
                svt=False
            if vphs[i] == vph and svt:
                line='  - '+filelists[i]+'\n'
                f.write(line)
                if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                    if ptt.exists(ob_path+'stds/'+filelists[i]+'.gz') == True:
                        call='gunzip -c '+ob_path+'stds/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                    else:
                        call='cp '+ob_path+'stds/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                    tools.sycall(call)
        f.close()

def create_obsLcbImage(data,vph,redux_path='',ob_path=''):
    """
    Creates a LcbImage observation file with the given name and path.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    nt=np.where((np.array(data['Type']) == 'IMAGE'))[0]
    if len(nt) == 0:
        print('No Object observations found')
        return
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    file=redux_path+'/obsresult-12VPH.yaml'.replace('VPH',vph)
    f=open(file,'w')
    f.write('id: 12VPH\n'.replace('VPH',vph))
    f.write('instrument: MEGARA\n')
    f.write('mode: MegaraLcbImage\n')
    f.write('images:\n')
    ct=0
    for i in range(0, len(filelists)):
        if vphs[i] == vph:
            line='  - '+filelists[i]+'\n'
            f.write(line)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == False:
                if ptt.exists(ob_path+'object/'+filelists[i]+'.gz') == True:
                    call='gunzip -c '+ob_path+'object/'+filelists[i]+'.gz > '+redux_path+'/data/'+filelists[i]
                else:
                    call='cp '+ob_path+'object/'+filelists[i]+' '+redux_path+'/data/'+filelists[i]
                tools.sycall(call)
            if ptt.exists(redux_path+'/data/'+filelists[i]) == True:    
                if ct == 0:
                    filetest=redux_path+'/data/'+filelists[i]
                    head=fits.getheader(filetest, 0)
                    insconfig=head['INSCONF']
                ct=ct+1
    f.close()
    return insconfig



def create_requirement(data,stds,vph,insconf,redux_path='',ob_path='',poly=False,calib_path='auxiliary/'):
    """
    Creates the create_requirement files with the given path.
    Parameters:
    data (dic): Dictionary from the observation file data
    stds (list): List of standard star names.
    vph (str): The Variable Passband Holography value.
    redux_path (str): The directory path where the file will be created.
    ob_path (str): The directory path where the observation files are located.
    """
    lampname=get_lamp(data,vph,calib_path=calib_path,redux_path=redux_path)
    tools.sycall('mkdir -p '+redux_path)
    tools.sycall('mkdir -p '+redux_path+'/data')
    for j in range(0, len(stds)):
        get_stdfile(stds[j],redux_path=redux_path,calib_path=calib_path)
        file=redux_path+'/requirements-'+stds[j]+'VPH.yaml'.replace('VPH',vph)
        f=open(file,'w')
        f.write('version: 1\n')
        f.write('products:\n')
        f.write('  MEGARA:\n')
        f.write('    INCONF: \n'.replace('INCONF',insconf))#ca3558e3-e50d-4bbc-86bd-da50a0998a48
        f.write("    - {id: 1, type: 'LinesCatalog', tags: {}, content: 'LAMP'}\n".replace('LAMP',lampname))
        f.write("    - {id: 2, type: 'MasterBias', tags: {}, content: 'master_bias.fits'}\n")
        f.write("    - {id: 3, type: 'TraceMap', tags: {}, content: 'master_tracesVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 4, type: 'MasterFiberFlat', tags: {}, content: 'master_fiberflatVPH.fits'}\n".replace('VPH',vph))
        f.write("    - {id: 5, type: 'WavelengthCalibration', tags: {}, content: 'master_wlcalibVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 6, type: 'MasterFiberFlatFrame', tags: {}, content: 'fiberflat_frameVPH.fits'}\n".replace('VPH',vph))
        f.write("    - {id: 7, type: 'ModelMap', tags: {}, content: 'master_modelVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 8, type: 'ReferenceSpectrumTable', tags: {}, content: 'mSTD.dat'}\n".replace('STD',stds[j]))
        f.write("    - {id: 9, type: 'ReferenceExtinctionTable', tags: {}, content: 'extintion_curve.txt'}\n")
        f.write("    - {id: 10, type: 'MasterSensitivity', tags: {}, content: 'master_sensitivityVPH_n.fits'}\n".replace('VPH',vph))
        if 'LR-U' in vph:
            f.write('requirements: \n')
            f.write('  MEGARA: \n')
            f.write('    default: \n')
            f.write('      default: \n')
            f.write('        MegaraArcImage: \n')
            f.write("          - {name: 'polynomial_degree', tags: {}, content: 5} \n")
            f.write("          - {name: 'nlines', tags: {}, content: [20, 20]} \n")
        else:
            f.write('requirements: {} \n')
        if poly:
            f.write('  MEGARA:\n')
            f.write('    MegaraArcCalibration:\n')
            f.write('      polynomial_degree: 2\n')
        f.close()
        file=redux_path+'/requirements-'+stds[j]+'VPHT.yaml'.replace('VPH',vph)
        f=open(file,'w')
        f.write('version: 1\n')
        f.write('products:\n')
        f.write('  MEGARA:\n')
        f.write('    INCONF: \n'.replace('INCONF',insconf))#ca3558e3-e50d-4bbc-86bd-da50a0998a48
        f.write("    - {id: 1, type: 'LinesCatalog', tags: {}, content: 'LAMP'}\n".replace('LAMP',lampname))
        f.write("    - {id: 2, type: 'MasterBias', tags: {}, content: 'master_bias.fits'}\n")
        f.write("    - {id: 3, type: 'TraceMap', tags: {}, content: 'master_tracesVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 4, type: 'MasterFiberFlat', tags: {}, content: 'master_fiberflatVPH.fits'}\n".replace('VPH',vph))
        f.write("    - {id: 5, type: 'WavelengthCalibration', tags: {}, content: 'master_wlcalibVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 6, type: 'MasterFiberFlatFrame', tags: {}, content: 'fiberflat_frameVPH.fits'}\n".replace('VPH',vph))
        f.write("    - {id: 7, type: 'ModelMap', tags: {}, content: 'master_modelVPH.json'}\n".replace('VPH',vph))
        f.write("    - {id: 8, type: 'ReferenceSpectrumTable', tags: {}, content: 'mSTD.dat'}\n".replace('STD',stds[j]))
        f.write("    - {id: 9, type: 'ReferenceExtinctionTable', tags: {}, content: 'extintion_curve.txt'}\n")
        if 'LR-U' in vph:
            f.write('requirements: \n')
            f.write('  MEGARA: \n')
            f.write('    default: \n')
            f.write('      default: \n')
            f.write('        MegaraArcImage: \n')
            f.write("          - {name: 'polynomial_degree', tags: {}, content: 5} \n")
            f.write("          - {name: 'nlines', tags: {}, content: [20, 20]} \n")
        else:
            f.write('requirements: {} \n')
        if poly:
            f.write('  MEGARA:\n')
            f.write('    MegaraArcCalibration:\n')
            f.write('      polynomial_degree: 2\n')
        f.close()
    if ptt.exists(redux_path+'/data/extintion_curve.txt') == False:
        call='cp '+calib_path+'extintion_curve.txt '+redux_path+'/data/extintion_curve.txt'
        tools.sycall(call)

def get_lamp(data,vph,calib_path='auxiliary/',redux_path=''):
    """
    Returns the lamp file based on the VPH value.
    Parameters:
    data (dic): Dictionary from the observation file data
    vph (str): The Variable Passband Holography value.

    Returns:
    str: The lamp file.
    """
    nt=np.where((np.array(data['Type']) == 'Arcs'))[0]
    if len(nt) == 0:
        print('No Arcs observations found')
        return
    root_p='megaradrp-calibrations-2018.1/LinesCatalog'
    dat=np.array(data['IMAGE'])
    dat1=np.array(data['VPH'])
    dat2=np.array(data['OSFILT'])
    dat4=np.array(data['ThAr1on'])
    #dat4=np.array(data['ThAr2on'])
    dat3=np.array(data['ThNe1on'])
    #dat6=np.array(data['ThAr3on'])
    #dat7=np.array(data['ThAr4on'])
    #dat8=np.array(data['ThAr5on'])
    filelists=dat[nt]
    vphs=dat1[nt]
    filt=dat2[nt]
    ThNe=dat3[nt]
    ThAr=dat4[nt]
    filelists,uin=np.unique(filelists,return_index=True)
    vphs=vphs[uin]
    filt=filt[uin]
    ThNe=ThNe[uin]
    ThAr=ThAr[uin]
    for i in range(0, len(filelists)):
        svt=False
        if filt[i] == 'BLUE' and ThAr[i] == 1:
            svt=True
        elif filt[i] == 'RED' and ThNe[i] == 1:
            svt=True
        else:
            svt=False
        if vphs[i] == vph and svt:
            if ThNe[i] == 1:
                pathL='ThNe'
            else:
                pathL='ThAr'
            file_name=vph+'_'+pathL+'.lis'
            pathT=calib_path+root_p+'/'+pathL+'/'+vph+'/'
            if ptt.exists(redux_path+'/data/'+file_name) == False:
                call='cp '+pathT+file_name+' '+redux_path+'/data/'+file_name
                tools.sycall(call)
    return file_name

def get_stdfile(std,redux_path='',calib_path='auxiliary/'):
    """
    Returns the standard star file based on the standard name.
    Parameters:
    std (str): The standard star name.

    Returns:
    str: The standard star file.
    """
    stdt=std.replace('HR','hr').replace('Feige','feige')
    root_p='esostandars/'
    if ptt.exists(calib_path+root_p+'ctiostan/mSTD.dat'.replace('STD',stdt)) == True:
        filestd1=calib_path+root_p+'ctiostan/mSTD.dat'.replace('STD',stdt)
        filestd2=calib_path+root_p+'ctiostan/fSTD.dat'.replace('STD',stdt)
    elif ptt.exists(calib_path+root_p+'okestan/mSTD.dat'.replace('STD',stdt)) == True:
        filestd1=calib_path+root_p+'okestan/mSTD.dat'.replace('STD',stdt)
        filestd2=calib_path+root_p+'okestan/fSTD.dat'.replace('STD',stdt)
    elif ptt.exists(calib_path+root_p+'hststan/mSTD.dat'.replace('STD',stdt)) == True:
        filestd1=calib_path+root_p+'hststan/mSTD.dat'.replace('STD',stdt)
        filestd2=calib_path+root_p+'hststan/fSTD.dat'.replace('STD',stdt)
    else:
        print('No standard star file found for: ',std)
        return
    if ptt.exists(redux_path+'/data/mSTD.dat'.replace('STD',std)) == False:
        call1='cp '+filestd1+' '+redux_path+'/data/mSTD.dat'.replace('STD',std)
        call2='cp '+filestd2+' '+redux_path+'/data/STD.dat'.replace('STD',std)
        tools.sycall(call1)
        tools.sycall(call2)
    return

def sort_megfiles(run,ob,obser_path='',path_redux='redux',calib_path='auxiliary/',forceAr=False):
    """
    Sorts the Megara files based on the observation file and creates necessary files.
    Parameters:
    run (str): The observation run identifier.
    ob (str): The observation identifier.
    obser_path (str): The directory path where the observation run are located.
    path_redux (str): The directory path where the reduced files will be created.
    calib_path (str): The directory path where the calibration files are located.
    """
    ob_path=obser_path+'OBob'.replace('ob',ob)+'/'
    rootf='RUN_ob_qc.txt'.replace('ob',ob).replace('RUN',run)
    print('Reading Megara observation file: ', rootf)
    data,hdr=read_obs(rootf,path=ob_path)
    obj_list=get_obj(data)
    vph_list=get_vph(data)
    std_list=get_std(data)
    insconfig=[]
    vphext,exp_list=get_exptime(data)
    print('Creating MegaraBiasImage file')
    create_obsBIAS(data,redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraTraceMap files for VPH: ', vph_list[i])
        create_obsTraceMap(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraArcCalibration files for VPH: ', vph_list[i])
        create_obsArcCalibration(data,vph_list[i],redux_path=path_redux,ob_path=ob_path,forceAr=forceAr)
    for i in range(0, len(vph_list)):
        print('Creating MegaraFiberFlatImage files for VPH: ', vph_list[i])
        create_obsFiberFlatImage(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraModelMap files for VPH: ', vph_list[i])
        create_obsModelMap(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraLcbStdStar files for VPH: ', vph_list[i])
        create_obsLcbStdStar(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraLcbImageStd files for VPH: ', vph_list[i])
        create_obsLcbImageStd(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
    for i in range(0, len(vph_list)):
        print('Creating MegaraLcbImage files for VPH: ', vph_list[i])
        insconf=create_obsLcbImage(data,vph_list[i],redux_path=path_redux,ob_path=ob_path)
        insconfig.extend([insconf])
    for i in range(0, len(vph_list)):
        print('Creating Requirement files for VPH: ', vph_list[i])
        create_requirement(data,std_list,vph_list[i],insconfig[i],redux_path=path_redux,ob_path=ob_path,calib_path=calib_path)
    return obj_list,std_list,vph_list,vphext,exp_list



def fix_lru_wavelength_calibration(
    filename,
    backup=True,
    bad_fibers=(1, 2, 604)
):
    """
    Fix known invalid LR-U wavelength solutions in a MEGARADRP
    WavelengthCalibration JSON file.

    Fibers 1 and 2:
        Copy the solution from fiber 3.

    Fiber 604:
        Interpolate the solution using fibers 603 and 605.

    Parameters
    ----------
    filename : str
        Path to master_wlcalibLR-U.json.

    backup : bool, optional
        If True, create filename + '.bak' before modifying.

    bad_fibers : tuple, optional
        Fibers to mark in error_fitting.

    Returns
    -------
    dict
        Modified wavelength calibration dictionary.
    """

    if backup:
        shutil.copy2(filename, filename + ".bak")

    with open(filename, "r") as f:
        data = json.load(f)

    by_fib = {
        item["fibid"]: item
        for item in data["contents"]
    }

    # ---------------------------------------------------------
    # Fibers 1 and 2:
    # copy wavelength solution from fiber 3
    # ---------------------------------------------------------

    for fib in (1, 2):
        by_fib[fib]["solution"] = copy.deepcopy(
            by_fib[3]["solution"]
        )

    # ---------------------------------------------------------
    # Fiber 604:
    # interpolate between fibers 603 and 605
    # ---------------------------------------------------------

    s603 = by_fib[603]["solution"]
    s604 = by_fib[604]["solution"]
    s605 = by_fib[605]["solution"]

    # Polynomial coefficients
    s604["coeff"] = [
        0.5 * (a + b)
        for a, b in zip(
            s603["coeff"],
            s605["coeff"]
        )
    ]

    # Linear wavelength solution
    for key in (
        "crpix",
        "crval",
        "cdelt",
        "crmin",
        "crmax"
    ):
        s604["cr_linear"][key] = 0.5 * (
            s603["cr_linear"][key]
            + s605["cr_linear"][key]
        )

    # Diagnostic quantities
    s604["residual_std"] = 0.5 * (
        s603["residual_std"]
        + s605["residual_std"]
    )

    s604["npoints_eff"] = int(
        round(
            0.5 * (
                s603["npoints_eff"]
                + s605["npoints_eff"]
            )
        )
    )

    s604["features"] = []

    # ---------------------------------------------------------
    # Mark fibers as problematic/interpolated
    # ---------------------------------------------------------

    previous_errors = set(
        data.get("error_fitting", [])
    )

    data["error_fitting"] = sorted(
        previous_errors | set(bad_fibers)
    )

    # ---------------------------------------------------------
    # Validate result
    # ---------------------------------------------------------

    invalid = []

    for item in data["contents"]:
        fibid = item["fibid"]
        cr_linear = item["solution"]["cr_linear"]

        if cr_linear["cdelt"] <= 0:
            invalid.append(fibid)

    if invalid:
        raise ValueError(
            "Invalid wavelength solutions remain "
            "for fibers: {}".format(invalid)
        )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

    #return data    

import json
import numpy as np


def generate_lru_healing(
    traces_file,
    output_file="healing_LRU.yaml",
    offset=30.0,
    xpoints=(100, 500, 1000, 1500, 2000,
             2500, 3000, 3500, 4000),
    poldeg=5,
    refit=True
):

    with open(traces_file, "r") as f:
        data = json.load(f)

    traces = data["contents"]

    blocks = []
    skipped = []

    for item in traces:

        fibid = item["fibid"]

        coeff = np.asarray(
            item.get("fitparms", []),
            dtype=float
        )

        # Skip fibers without a valid trace polynomial
        if coeff.size == 0:
            print(
                "WARNING: empty fitparms for fiber {}".format(
                    fibid
                )
            )
            skipped.append(fibid)
            continue

        xstart = item.get("start", 4)
        xstop = item.get("stop", 4092)

        user_points = []

        for x in xpoints:

            if x < xstart or x > xstop:
                continue

            y = np.polynomial.polynomial.polyval(
                float(x),
                coeff
            )

            y += offset

            # Extra safety
            if not np.isfinite(y):
                print(
                    "WARNING: non-finite trace for fiber {} at x={}".format(
                        fibid, x
                    )
                )
                continue

            user_points.append(
                (float(x), float(y))
            )

        if len(user_points) < poldeg + 1:

            print(
                "WARNING: insufficient points for fiber {}: {}".format(
                    fibid,
                    len(user_points)
                )
            )

            skipped.append(fibid)
            continue

        block = []

        block.append(
            "description: fit_through_user_points"
        )
        block.append(
            "fibid: {}".format(fibid)
        )
        block.append(
            "poldeg: {}".format(poldeg)
        )
        block.append(
            "xstart: {}".format(xstart)
        )
        block.append(
            "xstop: {}".format(xstop)
        )
        block.append("user_points:")

        for x, y in user_points:
            block.append(
                "  - [{:.1f}, {:.6f}]".format(
                    x, y
                )
            )

        block.append(
            "refit: {}".format(
                "True" if refit else "False"
            )
        )

        blocks.append("\n".join(block))

    with open(output_file, "w") as f:

        if blocks:
            f.write("\n---\n".join(blocks))
            f.write("\n")

    print("")
    print(
        "Generated {} healing blocks in {}".format(
            len(blocks),
            output_file
        )
    )

    if skipped:
        print(
            "Skipped {} fibers: {}".format(
                len(skipped),
                skipped
            )
        )

#    return output_file