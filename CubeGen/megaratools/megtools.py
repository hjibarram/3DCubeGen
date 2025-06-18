import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
import CubeGen.tools.tools as tools

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

def megarafiber_pos(hdr):
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
        
    #import matplotlib.pyplot as plt
    #plt.plot(x_posf*psc,y_posf*psc,'o')
    #plt.show()
    x_ifu=x_posf*psc
    y_ifu=y_posf*psc
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
    x_ifu,y_ifu,fib_idt,fib_ids=megarafiber_pos(hdr1)
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


def gen_sensf(wav_i,res_i,wave,flux,res_f,at_ext,minwave,maxwave,hdr2,path_data='data',vph='B',scalefact=1.0):
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
    return flux,s