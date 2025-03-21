from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy import units as u
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
from tqdm.notebook import tqdm
import os.path as ptt
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
import CubeGen.tools.tools as tools

def coadd_cube(nameR,nameF,path='',id_l=['0','1'],error=False,nsplit=0,spt=[0,0]):
    n_slides=len(id_l)
    cube_list=[]
    hdr_list=[]
    size_list=[]
    if error:
        cubeE_list=[]    
    pbar=tqdm(total=n_slides)  
    for i in range(0, n_slides):
        cube_file=path+'/'+nameR.replace('id',id_l[i])
        [cube, hdr]=fits.getdata(cube_file, 0, header=True)
        if i == 0:
            nzo,nxo,nyo=cube.shape
            #print(nzo,nxo,nyo)
            if nsplit > 1:
                nx_list=[]
                ny_list=[]
                dxt=int(nxo/nsplit)
                dyt=int(nyo/nsplit)
                for kt in range(0, nsplit):
                    nx_list.extend([dxt*kt])
                    ny_list.extend([dyt*kt])
                nx_list.extend([nxo])
                ny_list.extend([nyo])
                nx2o=nx_list[spt[0]+1]
                ny2o=ny_list[spt[1]+1]
                nx1o=nx_list[spt[0]]
                ny1o=ny_list[spt[1]]
                labf='_p'+str(spt[0])+str(spt[1])
            else:
                nx2o=nxo
                ny2o=nyo
                nx1o=0
                ny1o=0
                labf=''
            #print(nx2o,nx1o,nxo)
        nz,nx,ny=cube.shape        
       # cube=cube[:,nx1o:nx2o,ny1o:ny2o]
        cubeT=np.copy(cube[:,nx1o:nx2o,ny1o:ny2o])
        del cube
        cube=np.copy(cubeT)
        del cubeT
        if error:
            cubeE=fits.getdata(cube_file, 1, header=False)
          #  cubeE=cubeE[:,nx1o:nx2o,ny1o:ny2o]
            cubeT=np.copy(cubeE[:,nx1o:nx2o,ny1o:ny2o])
            del cubeE
            cubeE=np.copy(cubeT)
            del cubeT
            cubeE_list.extend([cubeE])
        cube_list.extend([cube])
        hdr_list.extend([hdr])
        size_list.extend([[nz,nx,ny]])
        pbar.update(1)
    pbar.close()
    
    
    wcs_list=[]
    wave_list=[]
    dpix_list=[]
    cdelt_list=[]
    for i in range(0, n_slides):
        nz=size_list[i][0]
        nx=size_list[i][1]
        ny=size_list[i][2]
        wcs=WCS(hdr_list[i])
        wcs=wcs.celestial
        wcs_list.extend([wcs])
        crpix=hdr_list[i]["CRPIX3"]
        try:
            cdelt=hdr_list[i]["CD3_3"]
        except:
            cdelt=hdr_list[i]["CDELT3"]
        cdelt_list.extend([cdelt])    
        crval=hdr_list[i]["CRVAL3"]
        wave=crval+cdelt*(np.arange(nz)+1-crpix)      
        wave_list.extend([wave])
        try:
            dx=np.sqrt((hdr_list[i]['CD1_1'])**2.0+(hdr_list[i]['CD1_2'])**2.0)*3600.0
            dy=np.sqrt((hdr_list[i]['CD2_1'])**2.0+(hdr_list[i]['CD2_2'])**2.0)*3600.0
        except:
            try:
                dx=hdr_list[i]['CD1_1']*3600.0
                dy=hdr_list[i]['CD2_2']*3600.0
            except:
                dx=hdr_list[i]['CDELT1']*3600.
                dy=hdr_list[i]['CDELT2']*3600.
        dpix=(np.abs(dx)+np.abs(dy))/2.0    
        dpix_list.extend([dpix])
        if i == 0:
            cdelt0=cdelt
            dpix0=dpix
            min_wave=np.nanmin(wave)
            nz0=nz
            nx0=nx
            ny0=ny
        if i == n_slides-1:
            max_wave=np.nanmax(wave)
    
    n_pix=int(np.round((max_wave-min_wave)/cdelt0))
    crval=min_wave
    crpix=1
    waveF=crval+cdelt0*(np.arange(n_pix)+1-crpix)
    #IFU_coadd=np.zeros([n_pix,nx0,ny0])
    #IFU_coaddE=np.zeros([n_pix,nx0,ny0])
    #IFU_coaddB=np.zeros([n_pix,nx0,ny0],dtype=int)
    IFU_coadd=np.zeros([n_pix,nx2o-nx1o,ny2o-ny1o])
    IFU_coaddE=np.zeros([n_pix,nx2o-nx1o,ny2o-ny1o])
    IFU_coaddB=np.zeros([n_pix,nx2o-nx1o,ny2o-ny1o],dtype=int)
    
    cdelt_a=np.nanmax(np.array(cdelt_list))
    pbar=tqdm(total=(nx2o-nx1o))#(nx0))#  
    for i in range(nx1o, nx2o):#0, nx0):#
        for j in range(ny1o, ny2o):#0, ny0):#
            temp_spec=np.copy(IFU_coadd[:,i-nx1o,j-ny1o])
            if error:
                temp_specE=np.copy(IFU_coaddE[:,i-nx1o,j-ny1o]) 
            pix1=i
            pix2=j
            for k in range(0, n_slides-1):
                nzi=size_list[k][0]
                nxi=size_list[k][1]
                nyi=size_list[k][2]
                wcs0=wcs_list[k]
                if k == 0:
                    sky1=pixel_to_skycoord(pix1,pix2,wcs0)
                xpos0a,ypos0a=skycoord_to_pixel(sky1,wcs0)
                xpos0=int(np.round(xpos0a))
                ypos0=int(np.round(ypos0a)) 
                if xpos0 >= 0 and xpos0 <= nxi-1 and ypos0 >= 0 and ypos0 <= nyi-1:
                    cube=cube_list[k]
                    spec0=cube[:,xpos0-nx1o,ypos0-ny1o] 
                    del cube
                    if error:
                        cubeE=cubeE_list[k]
                        specE0=cubeE[:,xpos0-nx1o,ypos0-ny1o]
                        del cubeE
                    #if np.nansum(spec0) != 0:
                    #    spec0=cube_interpolB(spec0,xpos0a,ypos0a)
                    if cdelt_a > cdelt_list[k]:
                        spec0=tools.median_a(spec,lw=cdelt_a)
                        if error:
                            specE0=tools.median_a(specE,lw=cdelt_a)
                    #print(xpos0,ypos0,'POS0',i,j)
                else:
                    spec0=np.zeros(nzi)
                    if error:
                        specE0=np.ones(nzi)
                nz1=size_list[k+1][0]
                nx1=size_list[k+1][1]
                ny1=size_list[k+1][2]
                wcs1=wcs_list[k+1]
                xpos1a,ypos1a=skycoord_to_pixel(sky1,wcs1)
                xpos1=int(np.round(xpos1a))
                ypos1=int(np.round(ypos1a))
                if xpos1 >= 0 and xpos1 <= nx1-1 and ypos1 >= 0 and ypos1 <= ny1-1:
                    cube=cube_list[k+1]
                    spec1=cube[:,xpos1-nx1o,ypos1-ny1o]
                    del cube
                    if error:
                        cubeE=cubeE_list[k+1]
                        specE1=cubeE[:,xpos1-nx1o,ypos1-ny1o]
                        del cubeE
                    #if np.nansum(spec1) != 0:
                    #    spec1=cube_interpolB(spec1,xpos1a,ypos1a)
                    if cdelt_a > cdelt_list[k+1]:
                        spec1=tools.median_a(spec1,lw=cdelt_a)
                        if error:
                            specE1=tools.median_a(specE1,lw=cdelt_a)
                    #print(xpos1,ypos1,'POS1',i,j)
                else:
                    spec1=np.zeros(nz1)
                    if error:
                        specE1=np.ones(nz1)  
                waveB=wave_list[k]
                waveG=wave_list[k+1]
                nt1i=np.where(waveB <= np.nanmin(waveG))
                nt2s=np.where(waveG >= np.nanmax(waveB))
                ntFs=np.where((waveF >= np.nanmax(waveB)) & (waveF <= np.nanmax(waveG)))
                ntFi=np.where((waveF <= np.nanmin(waveG)) & (waveF >= np.nanmin(waveB)))
                if len(nt1i[0]) > 0:     
                    specBFi=interp1d(waveB[nt1i],spec0[nt1i],bounds_error=False,fill_value=0.)(waveF[ntFi])
                    temp_spec[ntFi]=specBFi
                    if error:
                        specEBFi=interp1d(waveB[nt1i],specE0[nt1i],bounds_error=False,fill_value=0.)(waveF[ntFi])
                        temp_specE[ntFi]=specEBFi
                if len(nt2s[0]) > 0:    
                    specGFs=interp1d(waveG[nt2s],spec1[nt2s],bounds_error=False,fill_value=0.)(waveF[ntFs])
                    temp_spec[ntFs]=specGFs
                    if error:
                        specEGFs=interp1d(waveG[nt2s],specE1[nt2s],bounds_error=False,fill_value=0.)(waveF[ntFs])
                        temp_specE[ntFs]=specEGFs
            IFU_coadd[:,i-nx1o,j-ny1o]=temp_spec
            if error:
                IFU_coaddE[:,i-nx1o,j-ny1o]=temp_specE 
        pbar.update(1)
    pbar.close()
            #nt_z=np.where(temp_spec != 0)
            #if len(nt_z[0]) > 0:
            #    temp_specM=conv(temp_spec,ke=5)
            #    temp_specE=np.abs(temp_spec-temp_specM)
            #    temp_specE=np.sqrt(conv(temp_specE**2.0,ke=50))
            #    IFU_coaddE[:,i,j]=temp_specE*0.1
            #    ntp=np.where(temp_spec == 0)
            #    if len(ntp[0]) > 0:
            #        IFU_coaddE[ntp,i,j]=0.002
            #else:
            #    IFU_coaddE[:,i,j]=1.0
            #nt_z=np.where(temp_spec == 0)    
            #if len(nt_z[0]) > 0:
            #    IFU_coaddB[nt_z,i,j]=1
                
    hdr0=hdr_list[0]            
    h1=fits.PrimaryHDU(IFU_coadd)
    h2=fits.ImageHDU(IFU_coaddE)
    h3=fits.ImageHDU(IFU_coaddB)
    if nsplit > 1:
        h4=fits.ImageHDU()
    h_k=h1.header
    keys=list(hdr0.keys())
    for i in range(0, len(keys)):
        h_k[keys[i]]=hdr0[keys[i]]
        h_k.comments[keys[i]]=hdr0.comments[keys[i]]
    h_k["CRPIX1"]=hdr0["CRPIX1"]-ny1o
    h_k["CRPIX2"]=hdr0["CRPIX2"]-nx1o
    h_k['CDELT3']=cdelt
    h_k['CRPIX3']=crpix
    h_k['CRVAL3']=crval
    h_k.update()
    h_t=h2.header
    for i in range(0, len(keys)):
        h_t[keys[i]]=hdr0[keys[i]]
        h_t.comments[keys[i]]=hdr0.comments[keys[i]]
    h_t['EXTNAME'] ='Error_cube'
    h_t["CRPIX1"]=hdr0["CRPIX1"]-ny1o
    h_t["CRPIX2"]=hdr0["CRPIX2"]-nx1o
    h_t['CDELT3']=cdelt
    h_t['CRPIX3']=crpix
    h_t['CRVAL3']=crval
    h_t.update()
    h_r=h3.header
    for i in range(0, len(keys)):
        h_r[keys[i]]=hdr0[keys[i]]
        h_r.comments[keys[i]]=hdr0.comments[keys[i]]
    h_r['EXTNAME'] ='BADPIXELMASK'
    h_r["CRPIX1"]=hdr0["CRPIX1"]-ny1o
    h_r["CRPIX2"]=hdr0["CRPIX2"]-nx1o
    h_r['CDELT3']=cdelt
    h_r['CRPIX3']=crpix
    h_r['CRVAL3']=crval
    h_r.update()    
    if nsplit > 1:
        h_k=h4.header
        keys=list(hdr.keys())
        for i in range(0, len(keys)):
            h_k[keys[i]]=hdr[keys[i]]
            h_k.comments[keys[i]]=hdr.comments[keys[i]]
        h_k['NX0']=nxo
        h_k['NY0']=nyo                                       
        h_k['NX1']=nx1o
        h_k['NY1']=ny1o                                       
        h_k['NX2']=nx2o
        h_k['NY2']=ny2o                                       
        h_k.update()
        hlist=fits.HDUList([h1,h2,h3,h4])
    else:
        hlist=fits.HDUList([h1,h2,h3])
    hlist.update_extend()
    nameF=nameF+labf
    hlist.writeto(path+'/'+nameF.replace('.fits.gz','')+'.fits', overwrite=True)
    tools.sycall('gzip -f '+path+'/'+nameF.replace('.fits.gz','')+'.fits') 
                  
    return