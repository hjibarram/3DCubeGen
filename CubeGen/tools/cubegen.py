from tqdm.notebook import tqdm
from tqdm import tqdm as tqdmT
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy import units as u
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
from astropy.time import Time
from scipy.interpolate import interp1d
import os
from multiprocessing import Pool
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
import numpy as np
import CubeGen.tools.tools as tools
import CubeGen.tools.kernel as kernel 
import os.path as ptt

def map_ifu(expnumL,nameF=None,notebook=True,use_slitmap=True,errors=True,cent=False,coord_ast=[0,0],coord_cen=[0,0],pbars=True,flu16=False,multiT=False,spec_range=(None,None),fac_sizeX=1.0,fac_sizeY=1.0,pix_s=18.5,sigm_s=18.5,alph_s=2.0,out_path='',agcam_dir='',redux_ver='1.0.2.dev0',redux_dir='',tilelist=['11111'],tileglist=['0011XX'],mjd=['0000'],scp=112.36748321030637,basename='lvmCFrame-NAME.fits',basenameC='lvmCube-NAME.fits',path_lvmcore=''):
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
        if errors:
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
        equinox=Time(2024.8, format='jyear')
        equinox_J2000 = Time('J2000')
        coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame='fk5', equinox=equinox_J2000, unit='deg')
        newcoord = coord.transform_to('icrs')
        new_ra_fib=newcoord.ra.deg*3600.0
        new_dec_fib=newcoord.dec.deg*3600.0
        ra_fib=new_ra_fib[nt]
        dec_fib=new_dec_fib[nt]
        Std_id=Std_id[nt]

        
        
        if use_slitmap == False:
            agcam_coadd = agcam_dir+'/'+mjd[i % len(mjd)]+'/coadds/'+'lvm.sci.coadd_s'+expnum+'.fits'
            if os.path.isfile(agcam_coadd):
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
                ra0t=coord_ast[0]
                dec0t=coord_ast[1]
                if ra0t == 0 and dec0t == 0:
                    ra0t=np.mean(ra_fib)/3600.0
                    dec0t=np.mean(dec_fib)/3600.0
                wt1 = WCS(naxis=2)    
                wt1.wcs.crpix = [100, 100]
                wt1.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
                wt1.wcs.crval = [ra0t,dec0t]
                wt1.wcs.ctype = ["RA---TAN", "DEC--TAN"]
                wt1.wcs.radesys = 'ICRS'
                #wt1.wcs.equinox = 'J2000'#2024.8    
            
            if use_slitmap == False:
                rac0=rac
                dec0=dec
            nx0,ny0=rss.shape
            wave0=crval+cdelt*(np.arange(ny0)+1-crpix)
            wave0=wave0/(1+helio/299792.458)
            crval0=crval/(1+helio/299792.458)
            wave_1,wave_2=spec_range
            if wave_1 and wave_2:#wave_1 > np.nanmin(wave0) and
                if wave_1 < np.nanmax(wave0) and wave_2 > wave_1:# and wave_2 < np.nanmax(wave0):
                    nt=np.where((wave0 >= wave_1) & (wave0 <= wave_2))[0]
                    wave0=wave0[nt]
                    rss=rss[:,nt]
                    if errors:
                        erss=erss[:,nt]
                    crval0=np.nanmin(wave0)
                    ny0=len(wave0)
                else:
                    print('The wave Range is not well defined')  
                    return
            nfib0=len(Std_id)  
            rss_f=np.zeros([nfib0*nlt,ny0])
            rss_f[0:nfib0,:]=rss[Std_id,:]
            if errors:
                rss_ef=np.zeros([nfib0*nlt,ny0])
                rss_ef[0:nfib0,:]=erss[Std_id,:]
            x_ifu_V=np.zeros([nfib0*nlt,ny0])
            y_ifu_V=np.zeros([nfib0*nlt,ny0])
            if use_slitmap:
                x_ifu_pix=np.zeros([nfib0*nlt,ny0])
                y_ifu_pix=np.zeros([nfib0*nlt,ny0])
            for k in range(0, ny0):
                if use_slitmap == False:
                    ra_fib, dec_fib=tools.make_radec(xp,yp,rac,dec,PA)
                x_ifu_V[0:nfib0,k]=ra_fib
                y_ifu_V[0:nfib0,k]=dec_fib
                if use_slitmap:
                    sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
                    x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
                    x_ifu_pix[0:nfib0,k]=x_pixel
                    y_ifu_pix[0:nfib0,k]=y_pixel
        else:
            nx,ny=rss.shape  
            wave=crval+cdelt*(np.arange(ny)+1-crpix)
            wave=wave/(1+helio/299792.458)
            for j in range(0, nfib0):
                rss_f[nfib0*i+j,:]=interp1d(wave,rss[Std_id[j],:],kind='linear',bounds_error=False)(wave0)
                if errors:
                    rss_ef[nfib0*i+j,:]=interp1d(wave,erss[Std_id[j],:],kind='linear',bounds_error=False)(wave0)
            for k in range(0, ny0):
                if use_slitmap == False:
                    ra_fib, dec_fib=make_radec(xp,yp,rac,dec,PA)
                x_ifu_V[nfib0*i:nfib0*(i+1),k]=ra_fib
                y_ifu_V[nfib0*i:nfib0*(i+1),k]=dec_fib  
                if use_slitmap:
                    sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
                    x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
                    x_ifu_pix[nfib0*i:nfib0*(i+1),k]=x_pixel
                    y_ifu_pix[nfib0*i:nfib0*(i+1),k]=y_pixel
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
            yot=(np.amax(y_ifu_V[:,0])+np.amin(y_ifu_V[:,0]))/2.0
            xot=(np.amax(x_ifu_V[:,0])+np.amin(x_ifu_V[:,0]))/2.0
        else:
            sky_coord = SkyCoord(ra=ra_cen, dec=dec_cen, frame="icrs", unit="deg")
            xot, yot = skycoord_to_pixel(sky_coord, wt1)
            xot=xot*pix_s
            yot=yot*pix_s
    else:
        yot=(np.amax(y_ifu_V[:,0])+np.amin(y_ifu_V[:,0]))/2.0
        xot=(np.amax(x_ifu_V[:,0])+np.amin(x_ifu_V[:,0]))/2.0
    skycor = pixel_to_skycoord(xot/pix_s,yot/pix_s,wt1)
    xat=skycor.ra.value
    yat=skycor.dec.value
    x_ifu_V=x_ifu_V-xot
    y_ifu_V=y_ifu_V-yot
    nw=len(wave0)
    ns=len(x_ifu_V[:,0])
    fibA=35.3
    thet=0.0

    if cent:
        nlx=int(round((np.amax(x_ifu_V[:,0])-np.amin(x_ifu_V[:,0])+1)/pix_s))
        nly=int(round((np.amax(y_ifu_V[:,0])-np.amin(y_ifu_V[:,0])+1)/pix_s))
    else:
        nlx=int(round((np.amax([np.amax(x_ifu_V[:,0]),-np.amin(x_ifu_V[:,0])])+1)*2/pix_s))
        nly=int(round((np.amax([np.amax(y_ifu_V[:,0]),-np.amin(y_ifu_V[:,0])])+1)*2/pix_s))
    nlx=int(nlx*fac_sizeX)
    nly=int(nly*fac_sizeY)
    if nlx== 0:
        nlx=1
    if nly== 0:
        nly=1

    wt = WCS(naxis=2)
    #wt.wcs.crpix = [nlx/2+1, nly/2+0]
    wt.wcs.crpix = [nlx-(nlx/2+(100-xot/pix_s))+1.5, (nly/2+(100-yot/pix_s))-0.5]
    wt.wcs.cdelt = np.array([-pix_s/3600.0, pix_s/3600.0])
    #wt.wcs.crval = [xat,yat]
    wt.wcs.crval = [ra0t,dec0t]
    wt.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wt.wcs.radesys = 'ICRS'
    #wt.wcs.pc = [[np.cos(thet*np.pi/180.0), np.sin(thet*np.pi/180.0)],[-np.sin(thet*np.pi/180.0),  np.cos(thet*np.pi/180.0)]]
    #wt.wcs.equinox = 'J2000'


    #sky_coord = SkyCoord(ra=(x_ifu_V+xot)/3600., dec=(y_ifu_V+yot)/3600., frame="icrs", unit="deg")
    #x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt)
    #print(x_pixel.shape, nlx)

    ifu=np.zeros([nw,nly,nlx])
    ifu_e=np.ones([nw,nly,nlx])
    ifu_1=np.ones([nw,nly,nlx])
    ifu_m=np.zeros([nw,nly,nlx])
    xo=-nlx/2*pix_s
    yo=-nly/2*pix_s
    xi=xo
    xf=xo
    facto=(pix_s)**2.0/(np.pi*(fibA/2.0)**2.0)
    spec_ifu=rss_f*facto
    if errors:
        specE_ifu=rss_ef*facto 
    
    if pbars:
        if notebook:
            pbar=tqdm(total=nlx)
        else:     
            pbar=tqdmT(total=nlx)
    int_spect=np.zeros(nw)
    for i in range(0, nlx):
        xi=xf
        xf=xf+pix_s
        Rsp=np.sqrt((x_ifu_V[:,0]-(xf+xi)/2.0)**2.0)
        if sigm_s > fibA*3.5*2:
            ntp=np.where(Rsp <= (sigm_s/2.0))[0]
        else:    
            ntp=np.where(Rsp <= (fibA*3.5*2/2.0))[0]
        if multiT:
            nproc=3#3#cpu_count()
            with ThreadPool(nproc) as pool:
                if errors:
                    args=[(spec_ifu[ntp,:],specE_ifu[ntp,:],x_ifu_V[ntp,:],y_ifu_V[ntp,:],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc,errors) for npros in range(0, nproc)]                    
                else:
                    args=[(spec_ifu[ntp,:],None,x_ifu_V[ntp,:],y_ifu_V[ntp,:],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc,errors) for npros in range(0, nproc)]
                result_l = pool.map(kernel.task_wrapper, args)
        else:
            nproc=1
            npros=0
            result_l=[]
            if errors:
                args=(spec_ifu[ntp,:],specE_ifu[ntp,:],x_ifu_V[ntp,:],y_ifu_V[ntp,:],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc,errors)
            else:
                args=(spec_ifu[ntp,:],None,x_ifu_V[ntp,:],y_ifu_V[ntp,:],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc,errors)
            result_l.extend([kernel.task_wrapper(args)])
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
                ifu[:,j,nlx-(i+1)]=result[0][ct]
                if errors:
                    ifu_e[:,j,nlx-(i+1)]=result[1][ct]
                ct=ct+1
        
        if pbars:
            pbar.update(1)
    if pbars:
        pbar.close()

    if ptt.exists(out_path) == False:
        os.system('mkdir -p '+out_path)
    
    if flu16:
        ifu=ifu/1e-16
        ifu_e=ifu_e/1e-16#*100
    
    new_header = wt.to_header()
    h1=fits.PrimaryHDU(ifu,header=new_header)
    #h1=fits.PrimaryHDU(ifu)
    h2=fits.ImageHDU(ifu_e)
    h3=fits.ImageHDU(ifu_1)
    h4=fits.ImageHDU(ifu_m)
    head_list=[h1,h2,h3,h4]

    dx=0#+300.0/16.0/pix_s
    dy=0#+300.0/16.0/pix_s
    
    h=h1.header
    #h.extend(wt.to_header())
    keys=list(hdr0.keys())
    #for i in range(0, len(keys)):
    #    if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]:
    #        h[keys[i]]=hdr0[keys[i]]
    #        h.comments[keys[i]]=hdr0.comments[keys[i]]
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nlx
    h["NAXIS2"]=nly
    ##h["NDITER"]=(len(files),'Number of dither observations')
    ##h["BUNIT"]= ('1E-16 erg/s/cm^2','Unit of pixel value ' )
    ##h["OBJECT"]=hdr_0[0]['OBJECT']
    ##h["CTYPE"] = ("RA---TAN", "DEC--TAN")
    #h["CRVAL1"]=xat#xot/3600.0
    #h["CD1_1"]=-np.cos(thet*np.pi/180.)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    #h["CD1_2"]=-np.sin(thet*np.pi/180.)*pix_s/3600.0#*np.cos(yot/3600.0*np.pi/180.)
    #h["CRPIX1"]=nlx/2+0.5+dx
    #h["CTYPE1"]='RA---TAN'
    #h["CRVAL2"]=yat#yot/3600.0
    #h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    #h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    #h["CRPIX2"]=nly/2+0.5+dy
    #h["CTYPE2"]='DEC--TAN'
    #h["CUNIT1"]='deg     '                                           
    #h["CUNIT2"]='deg     '
    h["CDELT3"]=cdelt
    h["CD3_3"]=cdelt
    h["CRPIX3"]=crpix
    h["CRVAL3"]=crval0
    h["CUNIT3"]=('Angstrom','Units of coordinate increment and value    ')    
    h["CTYPE3"]=('WAVE    ','Air wavelength (linear) ')
    h["RADESYS"]='FK5     '
    h["OBJSYS"]='ICRS    '
    h["EQUINOX"]=2000.00
    h["IFUCON"]=(str(int(ns))+' ','NFibers')
    if flu16:
        h["BUNIT"]='10^-16 erg/s/cm^2'
    else:
        h["BUNIT"]='erg/s/cm^2'
    h.update() 
    hlist=fits.HDUList(head_list)
    hlist.update_extend()
    
    if nameF:
        file=out_path+basenameC.replace('NAME',nameF)
    else:
        file=out_path+basenameC.replace('NAME',expn)
    out_fit=file
    hlist.writeto(out_fit,overwrite=True)
    tools.sycall('gzip -f '+out_fit)