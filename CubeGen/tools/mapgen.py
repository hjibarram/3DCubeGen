from tqdm.notebook import tqdm
from tqdm import tqdm as tqdmT
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy import units as u
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
import os
from multiprocessing import Pool
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
import numpy as np
import CubeGen.tools.tools as tools
import CubeGen.tools.kernel as kernel 
import os.path as ptt

def gen_map(expnumL,nameF='MapLVM',notebook=True,use_slitmap=True,cent=False,coord_cen=[0,0],pbars=True,fac_sizeX=1.1,fac_sizeY=1.1,multiT=False,pix_s=18.5,zt=0,ki=5,sigm_s=18.5,alph_s=2.0,out_path='',agcam_dir='',redux_dir='',tilelist=['11111'],tileglist=['0011XX'],mjd=['0000'],redux_ver='0.1.1.dev0/1111/',scp=112.36748321030637,basename='lvmCFrame-NAME.fits',basenameC='lvmMap-NAME_TRA.fits',path_lvmcore=''):
    try:
        nlt=len(expnumL)
    except:
        nlt=1
    if pbars:    
        if notebook:
            pbar=tqdm(total=nlt)
        else:
            pbar=tqdmT(total=nlt)    
    for i in range(0, nlt):
        if nlt == 1:
            expnum=expnumL[i]
        else:
            expnum=expnumL[i]
        expn=str(int(expnum)).zfill(8)
        file=redux_dir+'/'+redux_ver+'/'+tileglist[i % len(tileglist)]+'/'+tilelist[i % len(tilelist)]+'/'+mjd[i % len(mjd)]+'/'+basename.replace('NAME',expn)
        hdr1=fits.getheader(file,0)
        [rss, hdr0]=fits.getdata(file,'FLUX', header=True)
        [ivar_rss, hdre]=fits.getdata(file,'IVAR', header=True)
        erss=1/np.sqrt(ivar_rss)    
        crpix=hdr0["CRPIX1"]
        cdelt=hdr0["CDELT1"]
        crval=hdr0["CRVAL1"]
        expT=float(hdr1['EXPTIME'])
        helio=float(hdr1['HIERARCH WAVE HELIORV_SCI'])
        
        hdu_list = fits.open(file)
        table_hdu = hdu_list['SLITMAP']
        table_data = table_hdu.data
        xp=table_data.field('xpmm')*scp
        yp=table_data.field('ypmm')*scp
        Std_id=table_data.field('fiberid')-1
        ra_fib=table_data.field('ra')*3600.0
        dec_fib=table_data.field('dec')*3600.0
        typ=table_data.field('targettype')
        nt=np.where(typ == 'science')
        xp=xp[nt]
        yp=yp[nt]
        ra_fib=ra_fib[nt]
        dec_fib=dec_fib[nt]
        Std_id=Std_id[nt]

        if use_slitmap == False:
            agcam_coadd = agcam_dir+'/'+mjd[i % len(mjd)]+'/coadds/'+'lvm.sci.coadd_s'+expnum+'.fits'
            if True:#os.path.isfile(agcam_coadd):
                agcam_hdu = fits.open(agcam_coadd)
                agcam_hdr = agcam_hdu[1].header
                w = WCS(agcam_hdr)
                cen = w.pixel_to_world(2500,1000)
                rac = cen.ra.deg  #agcam_hdr['RAMEAS']
                dec = cen.dec.deg #agcam_hdr['DECMEAS']hdr0
                PA = agcam_hdr['PAMEAS'] - 180.
                agcam_hdu.close()
            else:
                rac=hdr1["POSCIRA"]
                dec=hdr1["POSCIDE"]
                PA=hdr1["POSCIPA"]
    
        if i == 0:
            if use_slitmap:
                wt1 = WCS(naxis=2)    
                wt1.wcs.crpix = [100, 100]
                wt1.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
                wt1.wcs.crval = [np.mean(ra_fib)/3600.0,np.mean(dec_fib)/3600.0]
                wt1.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            
            if use_slitmap == False:
                rac0=rac
                dec0=dec
            nx0,ny0=rss.shape
            wave0=crval+cdelt*(np.arange(ny0)+1-crpix)
            wave0=wave0/(1+helio/299792.458)
            #crval0=crval/(1+helio/299792.458)
            nfib0=len(Std_id)  
            rss_f=np.zeros([nfib0*nlt])
            rss_ef=np.zeros([nfib0*nlt])
            rss_fm=np.zeros([nfib0*nlt])
            rss_efm=np.zeros([nfib0*nlt])
            for j in range(0, nfib0):
                val1,val2,val3,nam=tools.band_spectra(wave0,rss[Std_id[j],:],k=ki,zt=zt)
                rss_f[j]=val2
                rss_fm[j]=val1
                val1,val2,val3,nam=tools.band_spectra(wave0,erss[Std_id[j],:],k=ki,zt=zt)
                rss_ef[j]=val2
                rss_efm[j]=val1
            x_ifu_V=np.zeros([nfib0*nlt])
            y_ifu_V=np.zeros([nfib0*nlt])
            if use_slitmap:
                x_ifu_pix=np.zeros([nfib0*nlt])
                y_ifu_pix=np.zeros([nfib0*nlt])
            if use_slitmap == False:
                ra_fib, dec_fib=tools.make_radec(xp,yp,rac,dec,PA)
            x_ifu_V[0:nfib0]=ra_fib
            y_ifu_V[0:nfib0]=dec_fib
            if use_slitmap:
                sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
                x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
                x_ifu_pix[0:nfib0]=x_pixel
                y_ifu_pix[0:nfib0]=y_pixel
        else:
            nx,ny=rss.shape  
            wave=crval+cdelt*(np.arange(ny)+1-crpix)
            wave=wave/(1+helio/299792.458)
            for j in range(0, nfib0):
                val1,val2,val3,nam=tools.band_spectra(wave,rss[Std_id[j],:],k=ki,zt=zt)
                rss_f[nfib0*i+j]=val2
                rss_fm[nfib0*i+j]=val1
                val1,val2,val3,nam=tools.band_spectra(wave,erss[Std_id[j],:],k=ki,zt=zt)
                rss_ef[nfib0*i+j]=val2
                rss_efm[nfib0*i+j]=val1
            if use_slitmap == False:    
                ra_fib, dec_fib=tools.make_radec(xp,yp,rac,dec,PA)
            x_ifu_V[nfib0*i:nfib0*(i+1)]=ra_fib#xp+rac*3600
            y_ifu_V[nfib0*i:nfib0*(i+1)]=dec_fib#yp+dec*3600  
            if use_slitmap:
                sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
                x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
                x_ifu_pix[nfib0*i:nfib0*(i+1)]=x_pixel
                y_ifu_pix[nfib0*i:nfib0*(i+1)]=y_pixel
        if pbars:
            pbar.update(1)     
    if pbars:
        pbar.close()
    if use_slitmap:
        y_ifu_V=y_ifu_pix*pix_s
        x_ifu_V=x_ifu_pix*pix_s
    if cent:
        ra_cen,dec_cen=coord_cen
        if ra_cen == 0:
            yot=(np.amax(y_ifu_V)+np.amin(y_ifu_V))/2.0
            xot=(np.amax(x_ifu_V)+np.amin(x_ifu_V))/2.0
        else:
            sky_coord = SkyCoord(ra=ra_cen, dec=dec_cen, frame="icrs", unit="deg")
            xot, yot = skycoord_to_pixel(sky_coord, wt1)
            xot=xot*pix_s
            yot=yot*pix_s
    else:    
        yot=(np.amax(y_ifu_V)+np.amin(y_ifu_V))/2.0
        xot=(np.amax(x_ifu_V)+np.amin(x_ifu_V))/2.0
    skycor = pixel_to_skycoord(xot/pix_s,yot/pix_s,wt1)
    xat=skycor.ra.value
    yat=skycor.dec.value
    x_ifu_V=x_ifu_V-xot
    y_ifu_V=y_ifu_V-yot
    nw=len(wave0)
    ns=len(x_ifu_V)
    fibA=35.3
    thet=0.0

    
    if cent:
        nlx=int(round((np.amax(x_ifu_V)-np.amin(x_ifu_V)+1)/pix_s))
        nly=int(round((np.amax(y_ifu_V)-np.amin(y_ifu_V)+1)/pix_s))
    else:
        nlx=int(round((np.amax([np.amax(x_ifu_V),-np.amin(x_ifu_V)])+1)*2/pix_s))
        nly=int(round((np.amax([np.amax(y_ifu_V),-np.amin(y_ifu_V)])+1)*2/pix_s))
    nlx=int(nlx*fac_sizeX)
    nly=int(nly*fac_sizeY)
    if nlx== 0:
        nlx=1
    if nly== 0:
        nly=1
    
    wt = WCS(naxis=2)
    wt.wcs.crpix = [nlx/2+1, nly/2+0]
    #wt.wcs.cdelt = np.array([-np.cos(thet*np.pi/180.0)*pix_s/3600.0*np.cos(yot/3600.0*np.pi/180.), np.cos(thet*np.pi/180.0)*pix_s/3600.0])
    wt.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
    wt.wcs.crval = [xat,yat]
    wt.wcs.ctype = ["RA---TAN", "DEC--TAN"]

     
    #sky_coord = SkyCoord(ra=x_ifu_V+xot, dec=y_ifu_V+yot, frame="icrs", unit="deg")
    #x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt)

    ifu=np.zeros([nly,nlx])
    ifu_e=np.ones([nly,nlx])
    ifuM=np.zeros([nly,nlx])
    ifuM_e=np.ones([nly,nlx])
    ifu_1=np.ones([nly,nlx])
    ifu_m=np.zeros([nly,nlx])
    xo=-nlx/2*pix_s
    yo=-nly/2*pix_s
    xi=xo
    xf=xo
    facto=(pix_s)**2.0/(np.pi*(fibA/2.0)**2.0)
    spec_ifu=rss_f*facto
    specE_ifu=rss_ef*facto 
    specM_ifu=rss_fm-2.5*np.log10(facto)+5.0*np.log10(pix_s)
    specEM_ifu=rss_efm-2.5*np.log10(facto)+5.0*np.log10(pix_s)
    if pbars:
        if notebook:
            pbar=tqdm(total=nlx)
        else:     
            pbar=tqdmT(total=nlx)
    for i in range(0, nlx):
        xi=xf
        xf=xf+pix_s
        Rsp=np.sqrt((x_ifu_V-(xf+xi)/2.0)**2.0)
        if sigm_s > fibA*3.5*2:
            ntp=np.where(Rsp <= (sigm_s/2.0))[0]
        else:    
            ntp=np.where(Rsp <= (fibA*3.5*2/2.0))[0]
        if multiT:
            nproc=3#3#cpu_count()
            with ThreadPool(nproc) as pool:
                args=[(spec_ifu[ntp],specE_ifu[ntp],specM_ifu[ntp],specEM_ifu[ntp],x_ifu_V[ntp],y_ifu_V[ntp],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc) for npros in range(0, nproc)]                    
                #args=[(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,i,npros,nproc,wt,xot,yot) for npros in range(0, nproc)]                    
                result_l = pool.map(kernel.task_wrappermap, args)
        else:
            nproc=1
            npros=0
            result_l=[]
            args=(spec_ifu[ntp],specE_ifu[ntp],specM_ifu[ntp],specEM_ifu[ntp],x_ifu_V[ntp],y_ifu_V[ntp],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc)
            #args=(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,i,npros,nproc,wt,xot,yot)
            result_l.extend([kernel.task_wrappermap(args)])
        for npros in range(0, nproc):
            result=result_l[npros]
            val=int(nly/nproc)
            a1=val*npros
            if npros < nproc-1:
                a2=val*(npros+1)
            else:
                a2=nly
            ct=0
            for j in range(a1, a2):
                ifu[j,nlx-(i+1)]=result[0][ct]
                ifu_e[j,nlx-(i+1)]=result[1][ct]
                ifuM[j,nlx-(i+1)]=result[2][ct]
                ifuM_e[j,nlx-(i+1)]=result[3][ct]
                ct=ct+1
        if pbars:        
            pbar.update(1)
    if pbars:
        pbar.close()

    if ptt.exists(out_path) == False:
        os.system('mkdir -p '+out_path)


    
    
    h1=fits.PrimaryHDU(ifu)
    h2=fits.ImageHDU(ifuM)
    h3=fits.ImageHDU(ifu_e)
    h4=fits.ImageHDU(ifuM_e)
    head_list=[h1,h2,h3,h4]

    dx=0
    dy=0
    h=h1.header
    keys=list(hdr0.keys())
    for i in range(0, len(keys)):
        if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]:
            h[keys[i]]=hdr0[keys[i]]
            h.comments[keys[i]]=hdr0.comments[keys[i]]
    h["EXTNAME"]='FLUX'
    h["NAXIS"]=2 
    h["NAXIS1"]=nx
    h["NAXIS2"]=ny
    h["CRVAL1"]=xat#xot/3600.0#hdr1['CRVAL1']
    h["CD1_1"]=-np.cos(thet*np.pi/180.0)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    h["CD1_2"]=-np.sin(thet*np.pi/180.0)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    h["CRPIX1"]=nlx/2+0.5+dx#nlx/2+dx
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=yat#yot/3600.0#hdr1['CRVAL2']
    h["CD2_1"]=-np.sin(thet*np.pi/180.0)*pix_s/3600.0
    h["CD2_2"]=np.cos(thet*np.pi/180.0)*pix_s/3600.0
    h["CRPIX2"]=nly/2+0.5+dy
    h["CTYPE2"]='DEC--TAN'
    h["CUNIT1"]='deg     '                                           
    h["CUNIT2"]='deg     '
    h["RADESYS"]='FK5     '
    h["OBJSYS"]='ICRS    '
    h["EQUINOX"]=2000.00
    h["IFUCON"]=(str(int(ns))+' ','NFibers')
    h["BUNIT"]='erg/s/cm^2'
    ht=h2.header
    for i in range(0, len(keys)):
        if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]:
            ht[keys[i]]=hdr0[keys[i]]
            ht.comments[keys[i]]=hdr0.comments[keys[i]]
    ht["EXTNAME"]='MAG'
    ht["NAXIS"]=2 
    ht["NAXIS1"]=nx
    ht["NAXIS2"]=ny
    ht["CRVAL1"]=xat#xot/3600.0
    h["CD1_1"]=-np.cos(thet*np.pi/180.0)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    h["CD1_2"]=-np.sin(thet*np.pi/180.0)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    ht["CRPIX1"]=nlx/2+0.5+dx
    ht["CTYPE1"]='RA---TAN'
    ht["CRVAL2"]=yat#yot/3600.0
    ht["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    ht["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    ht["CRPIX2"]=nly/2+0.5+dy
    ht["CTYPE2"]='DEC--TAN'
    ht["CUNIT1"]='deg     '                                           
    ht["CUNIT2"]='deg     '
    ht["RADESYS"]='FK5     '
    ht["OBJSYS"]='ICRS    '
    ht["EQUINOX"]=2000.00
    ht["IFUCON"]=(str(int(ns))+' ','NFibers')
    ht["BUNIT"]='ABmag/arcsec'
    hlist=fits.HDUList(head_list)
    hlist.update_extend()
    basenameC=basenameC.replace('TRA',nam)
    file=out_path+basenameC.replace('NAME',nameF)
    out_fit=file
    #print(out_fit)
    hlist.writeto(out_fit,overwrite=True)
    tools.sycall('gzip -f '+out_fit)