import numpy as np
from scipy.interpolate import interp1d
#from astropy import wcs
from astropy.io import fits
#import warnings
#warnings.filterwarnings("ignore")
#from astropy.coordinates import ICRS, Galactic, FK4, FK5
#from astropy.wcs.utils import skycoord_to_pixel
#from astropy.wcs.utils import pixel_to_skycoord

import CubeGen.tools.tools as tools
import CubeGen.megaratools.megtools as mtools
import CubeGen.megaratools.megkernel as mkernel

nameF=None,notebook=True,ofvel=[[0,0,0],[0,0,0],[0,0,0]],use_slitmap=True,errors=True,cent=False,coord_ast=[0,0],coord_cen=[0,0],pbars=True,flu16=False,multiT=False,spec_range=(None,None),fac_sizeX=1.0,fac_sizeY=1.0,pix_s=18.5,sigm_s=18.5,alph_s=2.0,out_path='',agcam_dir='',redux_ver='1.0.2.dev0',redux_dir='',tilelist=['11111'],tileglist=['0011XX'],mjd=['0000'],scp=112.36748321030637,basename='lvmCFrame-NAME.fits',basenameC='lvmCube-NAME.fits',path_lvmcore=''
def megmap_ifu(nameF=None):
    """
    Generate a cube from MEGARA IFU data.
    
    This function reads MEGARA IFU data, processes it, and generates a cube with the specified wavelength range.
    It handles multiple observations and applies ADR corrections.
    
    Returns:
        None
    """

#files=['obsid5','obsid12']
files=['obsid12R']
if len (files) > 1:
    thiter=True
else:
    thiter=False
plot_sky=False

#wave_inf=4971#6109#4971#6109#7250#6109
#wave_sup=5435#7199#5435#7199#8600#7199#7099
wave_inf=4350#6109
wave_sup=5200#5196#7199#7099
wave_inf=6117
wave_sup=7199
#wave_inf=5170#5#5143#7250#6109
#wave_sup=6164#8600#7199#7099
dir_o=''
data_0=[]
hdr_0=[]
hdr_1=[]
for ii in range(0, len(files)):
    file='../'+files[ii]+'_results/final_rss.fits'
    [data, hdr]=fits.getdata(file,0, header=True)
    [data1, hdr1]=fits.getdata(file,1, header=True)
    print('Processing '+hdr['OBJECT'])
    n_fib,n_wave=data.shape
    print(data.shape)
    if plot_sky == True:
        flux=0
        for i in range(0, len(fib_ids)):
            flux=flux+data[np.int(fib_ids[i])-1,:]
        flux=flux/np.float(len(fib_ids))
        import matplotlib.pyplot as plt
        plt.xlim(wave_inf,wave_sup)
        plt.plot(wave1,flux)
        plt.show()
    if ii == 0:
        x_ifu,y_ifu,fib_idt,fib_ids=mtools.megarafiber_pos(hdr1)
        crval=hdr['CRVAL1']
        cdelt=hdr['CDELT1']
        crpix=hdr['CRPIX1']
        vel=mtools.get_radvel(hdr)
        crval=crval/(1+vel)
        cdelt=cdelt/(1+vel)
        wave1=crval+np.arange(n_wave)*cdelt
        nt=np.where((wave1 >= wave_inf) & (wave1 <= wave_sup))
        wave=wave1[nt]
        crval=wave[0]
        spec_ifu_t=data[:,nt[0]]/1e-16#*1.4575#*1.204
        outf=hdr['OBJECT']+'_R'
        spec_ifu=np.zeros([len(x_ifu),len(wave)])
        x_ifu_V=np.zeros([len(x_ifu),len(wave)])
        y_ifu_V=np.zeros([len(x_ifu),len(wave)])
        R2,R=mtools.get_adr(hdr,wave)
        Rt=np.zeros([2,len(wave)])
        Rt[0,:]=0
        Rt[1,:]=R#*10
        R_adr=np.dot(R2,Rt)
        print(Rt[1,0]*10,Rt[1,len(wave)-1])
        #print wave
        for i in range(0, len(x_ifu)):
            fib=int(fib_idt[i])-1
            spec_ifu[i,:]=spec_ifu_t[fib,:]
            x_ifu_V[i,:]=-R_adr[0,:]+x_ifu[i]
            y_ifu_V[i,:]=-R_adr[1,:]+y_ifu[i]
        data=data[:,nt[0]]
        hdr['CRVAL1']=crval
    else:
        x_ifu_d,y_ifu_d,fib_idt,fib_ids=mtools.megarafiber_pos(hdr1)
        wave1=hdr['CRVAL1']+np.arange(n_wave)*hdr['CDELT1']
        spec_ifu_d=np.zeros([len(x_ifu_d),len(wave)])
        
        x_ifu_dV=np.zeros([len(x_ifu_d),len(wave)])
        y_ifu_dV=np.zeros([len(x_ifu_d),len(wave)])
        R2,R=get_adr(hdr,wave)
        Rt=np.zeros([2,len(wave)])
        Rt[0,:]=0
        Rt[1,:]=R
        R_adr=np.dot(R2,Rt)
        
        for i in range(0, len(x_ifu_d)):
            fib=np.int(fib_idt[i])-1
            spec_ifu_d[i,:]=interp1d(wave1,data[fib,:],bounds_error=False,fill_value=0.)(wave)
            
            x_ifu_dV[i,:]=-R_adr[0,:]+x_ifu_d[i]
            y_ifu_dV[i,:]=-R_adr[1,:]+y_ifu_d[i]
            
        spec_ifu_d=spec_ifu_d/1e-16
        x_ifu=np.append(x_ifu,x_ifu_d)
        y_ifu=np.append(y_ifu,y_ifu_d)
        spec_ifu=np.append(spec_ifu,spec_ifu_d,axis=0)
        
        x_ifu_V=np.append(x_ifu_V,x_ifu_dV,axis=0)
        y_ifu_V=np.append(y_ifu_V,y_ifu_dV,axis=0)
        
        nt=np.where((wave1 >= wave_inf) & (wave1 <= wave_sup))
        data=data[:,nt[0]]
        hdr['CRVAL1']=wave1[nt][0]
    print(spec_ifu.shape,ii)
    data_0.extend([data])
    hdr_0.extend([hdr])
    hdr_1.extend([hdr1])




yot=(np.amax(y_ifu)+np.amin(y_ifu))/2.0
xot=(np.amax(x_ifu)+np.amin(x_ifu))/2.0
y_ifu=y_ifu-yot
x_ifu=x_ifu-xot
x_ifu_V=x_ifu_V-xot
y_ifu_V=y_ifu_V-yot
print(xot/3600.0,yot/3600.0)
print(y_ifu.shape)
nw=len(wave)
ns=len(x_ifu)
pix_s=0.35
fibA=0.000336666666666667/2*3600.0#*0.82#0.62#4.2
thet=0.0
nl=int(round((np.amax([np.amax(x_ifu),-np.amin(x_ifu),np.amax(y_ifu),-np.amin(y_ifu)])+1)*2/pix_s))

ifu=np.zeros([nw,nl,nl])
ifuE=np.ones([nw,nl,nl])
ifu_1=np.ones([nw,nl,nl])
ifu_m=np.zeros([nw,nl,nl])
xo=-nl/2*pix_s
yo=-nl/2*pix_s
xi=xo
xf=xo
facto=(pix_s)**2.0/(np.pi*(fibA/2.0)**2.0)
spec_ifu=spec_ifu*facto
sigm_s=pix_s
alph_s=2.0

for i in range(0, nl):
    xi=xf
    xf=xf+pix_s
    yi=yo
    yf=yo
    ifu,ifuE=mkernel.kernel_int(ifu,ifuE,spec_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yi,yf,xi,xf,nw,nl,i,erroF=False)

h1=fits.PrimaryHDU(ifu)
h2=fits.ImageHDU(ifuE)
h3=fits.ImageHDU(ifu_1)
h4=fits.ImageHDU(ifu_m)
head_list=[h1,h2,h3,h4]

for ii in range(0, len(files)):
    ifu_1=data_0[ii]
    ifu_m=[0]
    hdr_t=hdr_0[ii]
    hdr1_t=hdr_1[ii]
    h3=fits.ImageHDU(ifu_1)
    h4=fits.ImageHDU(ifu_m)
    ht=h3.header
    keys=list(hdr_t.keys())
    for i in range(0, len(keys)):
        if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]: 
            ht[keys[i]]=hdr_t[keys[i]]
            ht.comments[keys[i]]=hdr_t.comments[keys[i]]
    ht.update()
    
    ht1=h4.header
    keys=list(hdr1_t.keys())
    for i in range(0, len(keys)):
        if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]: 
            ht1[keys[i]]=hdr1_t[keys[i]]
            ht1.comments[keys[i]]=hdr1_t.comments[keys[i]]
    ht1['EXTNAME'] ='FIBERS_'+str(ii) 
    ht1.update()
    head_list.extend([h3])
    head_list.extend([h4])

dx=0
dy=0

h=h1.header

keys=list(hdr.keys())
for i in range(0, len(keys)):
    if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]:
        h[keys[i]]=hdr[keys[i]]
        h.comments[keys[i]]=hdr.comments[keys[i]]
del h["CDELT1"]
del h["CDELT2"]

h["NAXIS"]=3
h["NAXIS3"]=nw 
h["NAXIS1"]=nl
h["NAXIS2"]=nl
h["NDITER"]=(len(files),'Number of dither observations')
h["BUNIT"]= ('1E-16 erg/s/cm^2','Unit of pixel value ' )
h["OBJECT"]=hdr_0[0]['OBJECT']
h["CRVAL1"]=xot/3600.0#hdr1['CRVAL1']
h["CD1_1"]=-np.cos(thet*np.pi/180.)*pix_s/3600.
h["CD1_2"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
h["CRPIX1"]=nl/2+dx
h["CTYPE1"]='RA---TAN'
h["CRVAL2"]=yot/3600.0#hdr1['CRVAL2']
h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
h["CRPIX2"]=nl/2+0.5+dy
h["CTYPE2"]='DEC--TAN'
h['CUNIT1']='deg     '                                           
h['CUNIT2']='deg     '
h['CDELT3']=cdelt
h['CRPIX3']=crpix
h['CRVAL3']=crval
h['CUNIT3']=('Angstrom','Units of coordinate increment and value    ')    
h['CTYPE3']=('AWAV    ','Air wavelength (linear) ')
h['RADECSYS']='ICRS    '
h['SYSTEM']='FK5     '
h['EQUINOX']=2000.00
h['IFUCON']=(str(int(ns))+' ','NFibers')
h['UNITS']='1E-16 erg/s/cm^2'
h.update()
hlist=fits.HDUList(head_list)
hlist.update_extend()
out_fit=dir_o+outf+'.fits'
hlist.writeto(out_fit,overwrite=True)
dir_o1=dir_o.replace(" ",r'\ ')
out_fit1=dir_o1+outf+'.fits'
tools.sycall('gzip -f '+out_fit1)

