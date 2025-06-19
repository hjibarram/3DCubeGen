import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
import CubeGen.tools.tools as tools
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy import units as u

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
        #print 'FIB'+str(i+1)+'_X'
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
        x_ifu=x_posf*psc+hdr['CRVAL1']*3600.0
        y_ifu=y_posf*psc+hdr['CRVAL2']*3600.0 
    else:
    	x_ifu=x_posf*psc
        y_ifu=y_posf*psc
    if verbose:
        import matplotlib.pyplot as plt
        plt.plot(x_ifu/3600.,y_ifu/3600.,'o')
        plt.show() 
        #print(hdr['CRVAL2'])
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

def extract_cube(wave0,s,r,xo=0,yo=0,stdar_t='Feige32',path_ifu='',vph='B'):
    file=path_ifu+'/SPSTD_'+stdar_t+'_'+vph+'.fits.gz'
    [flux, hdrT]=fits.getdata(file, 0, header=True)
    nz,nx,ny=flux.shape
    rad=np.zeros([nx,ny])
    for i in range(0, nx):
        for j in range(0,ny):
            rt=np.sqrt((i-xo)**2.0+(j-yo)**2.0)*0.35
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

def get_rssflux_sens(r,xo=0,yo=0,path_block9='',path_sensfits='',path_data='data',vph='B',phase=0):
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
    xs=xo*0.35-7 #-npix/2*dpix IFU
    ys=yo*0.35-7 #-npix/2*dpix IFU
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
    #Second line
    #la1=6860#6884
    #la2=6884#7120
    #nt1=np.where((wav_i <= la1))
    #nt2=np.where((wav_i >= la2))
    #flutt1=res_i[nt1]
    #flutt2=res_i[nt2]
    #wavt1=wav_i[nt1]
    #wavt2=wav_i[nt2]
    #nt0=np.where((wave > la1) & (wave < la2))
    #flutt0=flux1[nt0]
    #wavt0=wave[nt0]
    #res_i=np.append(np.append(flutt1,flutt0),flutt2)
    #wav_i=np.append(np.append(wavt1,wavt0),wavt2)
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
    #nt3=np.where((wave_t >= 6860) & (wave_t <= 7400))#6884-7400
    #nt3=np.where((wave_t >= 8100) & (wave_t <= 8400))#6884-7400
    nt3=np.where((wave_t >= 6860) & (wave_t <= 7400))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    #nt3=np.where((wave_t >= 7617) & (wave_t <= 7720))#6884-7400
    nt3=np.where((wave_t >= 6860) & (wave_t <= 7400))#6884-7400
    resp_t1[nt3]=resp_t0[nt3]
    resp[nt]=resp_t1#model
    s=resp
    fits.writeto(path_data+'/master_sensitivity'+vph+'_n.fits', resp*scalefact, hdr2,overwrite=True)#*0.2*0.555#*0.3336#*0.864
    flux=flux/s
    return flux,s,res_f

def meg_spectra_outs(wave,flux,res_f,minwave,maxwave,vph='B'):
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
    plt.savefig('spec'+vph+'.jpg')#,dpi=1000)
    plt.close()
    return factor

def calib_spec(phase=0,scfact=1.0,xo=0,yo=0,vph='B',stdar_t='Feige32'):
    r=4.0#Aperture radius
    path_data='data'
    path_block9='obsid9'+vph+'_results/'
    path_sensfits='obsid9'+vph+'_results/'
    if phase == 0:
        cube=False
        gen_hr=False
        gen_sensf=True
        scalefact=1.0
        fergs=True
        stdT=''
    if phase == 1:
        cube=False
        gen_hr=True
        gen_sensf=False
        scalefact=1.0
        fergs=True
        stdT=''
    if phase == 2:
        cube=False
        gen_hr=False
        gen_sensf=True
        scalefact=1.0
        fergs=False
        stdT='_hr'+vph
    if phase == 3 or phase == 5:
        cube=True
        gen_hr=False
        gen_sensf=False
        scalefact=1.0
        fergs=False
        stdT='_hr'+vph
    if phase == 4:
        cube=False
        gen_hr=False
        gen_sensf=True
        scalefact=scfact
        fergs=False
        stdT='_hr'+vph

    wav_i,res_i=read_standar(path_data=path_data,stdar_t=stdar_t,stdT=stdT,fergs=fergs)
    wave,flux,Xa,s,hdr2=get_rssflux_sens(r,xo=xo,yo=yo,path_block9=path_block9,path_sensfits=path_sensfits,path_data=path_data,vph=vph,phase=phase)
    if cube == True:
        flux,wave,s=extract_cube(wave,s,r,xo=xo,yo=yo,stdar_t=stdar_t,path_ifu='ifu',vph=vph)
    maxwave=np.round(np.nanmax(wave)-70)
    minwave=np.round(np.nanmin(wave)+70)
    Kvl=extintion_c(wave)
    at_ext=10.0**(-0.4*Xa*Kvl)
    flux1=flux/s/at_ext
    if gen_hr == True:
        wav_i,res_i=high_res_std(wav_i,res_i,wave,flux,flux1,path_data=path_data,stdar_t=stdar_t,vph=vph)
    res_f=interp1d(wav_i,res_i,bounds_error=False,fill_value=0.)(wave)
    if gen_sensf == True:
        flux,s,res_f=gen_sensf(wav_i,wave,flux,res_f,at_ext,minwave,maxwave,hdr2,path_data=path_data,vph=vph,scalefact=scalefact)
    factor=meg_spectra_outs(wave,flux,res_f,minwave,maxwave,vph=vph)
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

def get_adr(hdr,wave,lo=5000.0):
    mjd=hdr["MJD-OBS"]
    lng=hdr["LONGITUD"].replace("+","-")
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

def get_radvel(hdr):
    import astropy.units as u
    mjd=hdr["MJD-OBS"]
    lng=hdr["LONGITUD"].replace("+","-")
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