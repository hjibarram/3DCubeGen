from tqdm.notebook import tqdm
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
from scipy.spatial.distance import pdist, squareform

def gen_matrix(expnumL,multiT=False,errors=True,covana=False,nprocf=6,pix_s=18.5,fac_sizeX=1.1,fac_sizeY=1.1,ki=5,sigm_s=18.5,alph_s=2.0,verbose=True,agcam_dir='',redux_dir='',tilelist=['11111'],tileglist=['0011XX'],mjd=['0000'],redux_ver='1.1.1.dev0/',scp=112.36748321030637,basename='lvmCFrame-NAME.fits',path_lvmcore=''):
    try:
        nlt=len(expnumL)
    except:
        nlt=1
    if verbose:    
        pbar=tqdm(total=nlt)    
    for i in range(0, nlt):
        if nlt == 1:
            expnum=expnumL[i]
        else:
            expnum=expnumL[i]
        expn=str(int(expnum)).zfill(8)
        file=redux_dir+'/'+redux_ver+'/'+tileglist[i % len(tileglist)]+'/'+tilelist[i % len(tilelist)]+'/'+mjd[i % len(mjd)]+'/'+basename.replace('NAME',expn)
        #file=redux_dir+'/'+redux_ver+'/'+mjd[i % len(mjd)]+'/'+basename.replace('NAME',expn)
        hdr1=fits.getheader(file,0)
        [rss, hdr0]=fits.getdata(file,'FLUX', header=True)
        [ivar_rss, hdre]=fits.getdata(file,'IVAR', header=True)
        erss=1/np.sqrt(ivar_rss)    
        crpix=hdr0["CRPIX1"]
        cdelt=hdr0["CDELT1"]
        crval=hdr0["CRVAL1"]
        expT=float(hdr1['EXPTIME'])
        
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
    
        if i == 0:
            wt1 = WCS(naxis=2)    
            wt1.wcs.crpix = [100, 100]
            wt1.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
            wt1.wcs.crval = [np.mean(ra_fib)/3600.0,np.mean(dec_fib)/3600.0]
            wt1.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            
            nx0,ny0=rss.shape
            wave0=crval+cdelt*(np.arange(ny0)+1-crpix)
            wave0=wave0
            nfib0=len(Std_id)  
            rss_f=np.zeros([nfib0*nlt])
            rss_ef=np.zeros([nfib0*nlt])
            rss_fm=np.zeros([nfib0*nlt])
            rss_efm=np.zeros([nfib0*nlt])
            for j in range(0, nfib0):
                val1,val2,val3,nam=tools.band_spectra(wave0,rss[Std_id[j],:],k=ki)
                rss_f[j]=val2
                rss_fm[j]=val1
                val1,val2,val3,nam=tools.band_spectra(wave0,erss[Std_id[j],:],k=ki)
                rss_ef[j]=val2
                rss_efm[j]=val1
            x_ifu_V=np.zeros([nfib0*nlt])
            y_ifu_V=np.zeros([nfib0*nlt])
            x_ifu_pix=np.zeros([nfib0*nlt])
            y_ifu_pix=np.zeros([nfib0*nlt])
            x_ifu_V[0:nfib0]=ra_fib
            y_ifu_V[0:nfib0]=dec_fib
            sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
            x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
            x_ifu_pix[0:nfib0]=x_pixel
            y_ifu_pix[0:nfib0]=y_pixel
        else:
            nx,ny=rss.shape  
            wave=crval+cdelt*(np.arange(ny)+1-crpix)
            wave=wave
            for j in range(0, nfib0):
                val1,val2,val3,nam=tools.band_spectra(wave,rss[Std_id[j],:],k=ki)
                rss_f[nfib0*i+j]=val2
                rss_fm[nfib0*i+j]=val1
                val1,val2,val3,nam=tools.band_spectra(wave,erss[Std_id[j],:],k=ki)
                rss_ef[nfib0*i+j]=val2
                rss_efm[nfib0*i+j]=val1
            x_ifu_V[nfib0*i:nfib0*(i+1)]=ra_fib#xp+rac*3600
            y_ifu_V[nfib0*i:nfib0*(i+1)]=dec_fib#yp+dec*3600  
            sky_coord = SkyCoord(ra=ra_fib/3600.0, dec=dec_fib/3600.0, frame="icrs", unit="deg")
            x_pixel, y_pixel = skycoord_to_pixel(sky_coord, wt1)
            x_ifu_pix[nfib0*i:nfib0*(i+1)]=x_pixel
            y_ifu_pix[nfib0*i:nfib0*(i+1)]=y_pixel
        if verbose:
            pbar.update(1)     
    if verbose:
        pbar.close()
    y_ifu_V=y_ifu_pix*pix_s
    x_ifu_V=x_ifu_pix*pix_s
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


    nlx=int(round((np.amax([np.amax(x_ifu_V),-np.amin(x_ifu_V)])+1)*2/pix_s))
    nly=int(round((np.amax([np.amax(y_ifu_V),-np.amin(y_ifu_V)])+1)*2/pix_s))
    nlx=int(nlx*fac_sizeX)
    nly=int(nly*fac_sizeY)
    
    wt = WCS(naxis=2)
    wt.wcs.crpix = [nlx/2+1, nly/2+0]
    wt.wcs.cdelt = np.array([pix_s/3600.0, pix_s/3600.0])
    wt.wcs.crval = [xat,yat]
    wt.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    if covana:
        Wt=np.zeros([nly,nlx,ns])
        St=np.zeros([ns,ns])
        for i in range(0,ns):
            St[i,i]=rss_ef[i]
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
    
    
    if verbose:
        pbar=tqdm(total=nlx)
    for i in range(0, nlx):
        xi=xf
        xf=xf+pix_s
        Rsp=np.sqrt((x_ifu_V-(xf+xi)/2.0)**2.0)
        #ntp=np.where(Rsp <= (fibA*3.5*2/2.0))[0]
        if sigm_s > fibA*3.5*2:
            radiT=sigm_s/2.0
        else:    
            radiT=fibA*3.5*2/2.0
        ntp=np.where(Rsp <= (radiT))[0]    

        if covana:
            xtF=x_ifu_V[ntp]
            ytF=y_ifu_V[ntp]
            for j in range(0, nly):
                yi=yo+pix_s*j
                yf=yo+pix_s*(j+1)
                Wgt1=0
                Rspt=np.sqrt((ytF-(yf+yi)/2.0)**2.0)
                ntpt=np.where(Rspt <= radiT)[0]
                for k in range(0, len(xtF[ntpt])):
                    Rsp=np.sqrt((xtF[ntpt][k]-(xf+xi)/2.0)**2.0+(ytF[ntpt][k]-(yf+yi)/2.0)**2.0)
                    if Rsp <= (radiT): 
                        Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)
                        Wgt1=Wgt1+Wg
                if Wgt1 == 0:
                    Wgt1=1
            
                for k in range(0, len(xtF[ntpt])):
                    Rsp=np.sqrt((xtF[ntpt][k]-(xf+xi)/2.0)**2.0+(ytF[ntpt][k]-(yf+yi)/2.0)**2.0)
                    if Rsp <= (radiT): 
                        Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)/Wgt1
                        ifi=ntp[ntpt][k]
                        Wt[j,nlx-(i+1),ifi]=Wg

        if multiT:
            nproc=nprocf#3#cpu_count()
            with ThreadPool(nproc) as pool:
                args=[(spec_ifu[ntp],specE_ifu[ntp],specM_ifu[ntp],specEM_ifu[ntp],x_ifu_V[ntp],y_ifu_V[ntp],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc) for npros in range(0, nproc)]                    
                #args=[(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,i,npros,nproc,wt,xot,yot) for npros in range(0, nproc)]                    
                result_l = pool.map(kernel.task_wrappermatrix, args)
        else:
            nproc=1
            npros=0
            result_l=[]
            args=(spec_ifu[ntp],specE_ifu[ntp],specM_ifu[ntp],specEM_ifu[ntp],x_ifu_V[ntp],y_ifu_V[ntp],fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,npros,nproc)
            #args=(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nly,i,npros,nproc,wt,xot,yot)
            result_l.extend([kernel.task_wrappermatrix(args)])
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
        if verbose:        
            pbar.update(1)
    if verbose:
        pbar.close()
        
    #print(np.nansum(ifu))
    #sys.exit()
    h1=fits.PrimaryHDU(ifu)
    h2=fits.ImageHDU(ifu_e)
    dx=0.5
    dy=0#-1
    h=h1.header
    keys=list(hdr0.keys())
    for i in range(0, len(keys)):
        if not "COMMENT" in  keys[i] and not 'HISTORY' in keys[i]:
            h[keys[i]]=hdr0[keys[i]]
            h.comments[keys[i]]=hdr0.comments[keys[i]]
    #del h["CDELT1"]
    #del h["WCSAXES"]
    h["EXTNAME"]='FLUX'
    h["NAXIS"]=2 
    h["NAXIS1"]=nlx
    h["NAXIS2"]=nly
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

    


    if covana:
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
                    Wgt=Wg[:,a1:a2,:]     
                    args.extend([(St,Wgt,nly,ns,a1,a2)])                    
                result_l = pool.map(task_wrappercov1, args)
        else:
            nproc=1
            npros=0
            result_l=[]
            args=(St,Wg,nly,ns,0,nlx)
            result_l.extend([task_wrappercov1(args)])
        for npros in range(0, nproc):
            result=result_l[npros]
            val=int(nlx/nproc)
            a1=val*npros
            if npros < nproc-1:
                a2=val*(npros+1)
            else:
                a2=nlx
            out[:,a1:a2,:]=result[0]
        
        out2=np.zeros([nly,nlx,nly,nlx])
        distF=np.zeros([nly,nlx,nly,nlx])#nx
        #verbose=True
        indexT=np.array([(k,0) for k in range(nly)])
        Dq=pdist(indexT, metric='euclidean')
        #start = time()
        if verbose:
            pbar=tqdm(total=nlx)
        for i in range(0, nlx):
            Wgt=Wg[:,i,:]
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
                    result_l = pool.map(task_wrappercov2, args)
            else:
                nproc=1
                npros=0
                result_l=[]
                args=(St,Wgt,out,Dq,nly,i,0,nlx)
                result_l.extend([task_wrappercov2(args)])
            for npros in range(0, nproc):
                result=result_l[npros]
                val=int(nlx/nproc)
                a1=val*npros
                if npros < nproc-1:
                    a2=val*(npros+1)
                else:
                    a2=nlx
                out2[:,i,:,a1:a2]=result[0]
                distF[:,i,:,a1:a2]=result[1]
            if verbose:        
                pbar.update(1)
        if verbose:
            pbar.close()

    
    if covana:
        return h,ifu,ifu_e,ifuM,ifuM_e,Wt,St
    else:
        return h,ifu,ifu_e,ifuM,ifuM_e