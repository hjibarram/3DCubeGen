import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.wcs import WCS
import CubeGen.tools.tools as tools
import CubeGen.megaratools.megtools as mtools
import CubeGen.megaratools.megkernel as mkernel

def megmap_ifu(reduxL,nameF=None,errors=False,flu16=True,spec_range=(None,None),fac_sizeX=1.0,fac_sizeY=1.0,pix_s=0.35,sigm_s=0.35,alph_s=2.0,out_path='',redux_dir='',vph='R',scp=112.36748321030637,basename='final_rss.fits',basenameC='megCube-NAME.fits'):
    """
    Generate a cube from MEGARA IFU data.
    
    This function reads MEGARA IFU data, processes it, and generates a cube with the specified wavelength range.
    It handles multiple observations and applies ADR corrections.
    
    Returns:
        None
    """
    try:
        nlt=len(reduxL)
    except:
        nlt=1
    data_0=[]
    hdr_0=[]
    hdr_1=[]
    for ii in range(0, nlt):
        file=redux_dir+'/'+reduxL[ii]+'_results/'+basename
        [rss, hdr]=fits.getdata(file,0, header=True)
        [erss, hdr1]=fits.getdata(file,1, header=True)
        print('Processing '+hdr['OBJECT'])
        n_fib,ny0=rss.shape
        if ii == 0:
            outf=hdr['OBJECT']+'_R'
            x_ifu,y_ifu,fib_idt,fib_ids=mtools.megarafiber_pos(hdr1)
            crval=hdr['CRVAL1']
            cdelt=hdr['CDELT1']
            crpix=hdr['CRPIX1']
            vel=mtools.get_radvel(hdr)
            crval=crval/(1+vel)
            cdelt=cdelt/(1+vel)
            wave0=crval+np.arange(ny0)*cdelt
            wave_1,wave_2=spec_range
            if wave_1 and wave_2:
                if wave_1 < np.nanmax(wave0) and wave_2 > wave_1:
                    nt=np.where((wave0 >= wave_1) & (wave0 <= wave_2))[0]
                    wave0=wave0[nt]
                    rss=rss[:,nt]
                    if errors:
                        erss=erss[:,nt]
                    crval=np.nanmin(wave0)
                    ny0=len(wave0)
                else:
                    print('The wave Range is not well defined')  
                    return
            else:
                wave_2=np.round(np.nanmax(wave0)-70)
                wave_1=np.round(np.nanmin(wave0)+70)
                nt=np.where((wave0 >= wave_1) & (wave0 <= wave_2))[0]
                wave0=wave0[nt]
                rss=rss[:,nt]
                if errors:
                    erss=erss[:,nt]
                crval=np.nanmin(wave0)
                ny0=len(wave0)
            n_fib0=len(x_ifu)
            if errors:
                rss_ef=np.zeros([nfib0*nlt,ny0])
            rss_f=np.zeros([n_fib0*nlt,ny0])
            x_ifu_V=np.zeros([n_fib0*nlt,ny0])
            y_ifu_V=np.zeros([n_fib0*nlt,ny0])
            R2,R=mtools.get_adr(hdr,wave0)
            Rt=np.zeros([2,ny0])
            Rt[0,:]=0
            Rt[1,:]=R
            R_adr=np.dot(R2,Rt)
            for i in range(0, n_fib0):
                fib=int(fib_idt[i])-1
                rss_f[i,:]=rss[fib,:]
                if errors:
                    rss_ef[i,:]=erss[fib,:]
                x_ifu_V[i,:]=-R_adr[0,:]+x_ifu[i]
                y_ifu_V[i,:]=-R_adr[1,:]+y_ifu[i]
            hdr['CRVAL1']=crval
        else:
            x_ifu,y_ifu,fib_idt,fib_ids=mtools.megarafiber_pos(hdr1)
            wave=hdr['CRVAL1']+np.arange(ny0)*hdr['CDELT1']
            R2,R=get_adr(hdr,wave0)
            Rt=np.zeros([2,ny0])
            Rt[0,:]=0
            Rt[1,:]=R
            R_adr=np.dot(R2,Rt)
            for i in range(0, len(x_ifu)):
                fib=np.int(fib_idt[i])-1
                rss_f[nfib0*ii+i,:]=interp1d(wave,rss[fib,:],kind='linear',bounds_error=False)(wave0)
                if errors:
                    rss_ef[nfib0*ii+i,:]=interp1d(wave,erss[fib,:],kind='linear',bounds_error=False)(wave0)
                x_ifu_V[nfib0*ii+i,:]=-R_adr[0,:]+x_ifu[i]
                y_ifu_V[nfib0*ii+i,:]=-R_adr[1,:]+y_ifu[i]
            nt=np.where((wave1 >= wave_1) & (wave1 <= wave_2))
            rss=rss_f[nfib0*ii:nfib0*(ii+1),:]
            hdr['CRVAL1']=crval
        data_0.extend([rss])
        hdr_0.extend([hdr])
        hdr_1.extend([hdr1])
    yot=(np.amax(y_ifu_V[:,0])+np.amin(y_ifu_V[:,0]))/2.0
    xot=(np.amax(x_ifu_V[:,0])+np.amin(x_ifu_V[:,0]))/2.0
    x_ifu_V=x_ifu_V-xot
    y_ifu_V=y_ifu_V-yot
    nw=len(wave0)
    ns=len(x_ifu_V[:,0])
    pix_s=0.35
    fibA=0.000336666666666667/2*3600.0#
    thet=0.0
    nlx=int(round((np.amax([np.amax(x_ifu_V[:,0]),-np.amin(x_ifu_V[:,0])])+1)*2/pix_s))
    nly=int(round((np.amax([np.amax(y_ifu_V[:,0]),-np.amin(y_ifu_V[:,0])])+1)*2/pix_s))
    nlx=int(nlx*fac_sizeX)
    nly=int(nly*fac_sizeY)
    if nlx== 0:
        nlx=1
    if nly== 0:
        nly=1
    wt = WCS(naxis=2)
    wt.wcs.crpix = [nlx/2+0, nly/2+0.5]
    #wt.wcs.crpix = [nlx-(nlx/2+(100-xot/pix_s))+1.5, (nly/2+(100-yot/pix_s))-0.5]
    wt.wcs.cdelt = np.array([-pix_s/3600.0, pix_s/3600.0])
    #wt.wcs.crval = [xat,yat]
    wt.wcs.crval = [xot/3600.0,yot/3600.0]
    wt.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wt.wcs.radesys = 'ICRS'

    ifu=np.zeros([nw,nly,nlx])
    ifuE=np.ones([nw,nly,nlx])
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
    for i in range(0, nlx):
        xi=xf
        xf=xf+pix_s
        yi=yo
        yf=yo
        ifu,ifuE=mkernel.kernel_int(ifu,ifuE,spec_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yi,yf,xi,xf,nw,nly,i,erroF=False)
    if flu16:
        ifu=ifu/1e-16
        ifuE=ifuE/1e-16
    h1=fits.PrimaryHDU(ifu)
    h2=fits.ImageHDU(ifuE)
    h3=fits.ImageHDU(ifu_1)
    h4=fits.ImageHDU(ifu_m)
    head_list=[h1,h2,h3,h4]
    for ii in range(0, len(reduxL)):
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
    h["NAXIS1"]=nlx
    h["NAXIS2"]=nly
    h["NDITER"]=(len(reduxL),'Number of dither observations')
    h["BUNIT"]= ('1E-16 erg/s/cm^2','Unit of pixel value ' )
    h["OBJECT"]=hdr_0[0]['OBJECT']
    h["CRVAL1"]=xot/3600.0#hdr1['CRVAL1']
    h["CD1_1"]=-np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nlx/2+dx
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=yot/3600.0#hdr1['CRVAL2']
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nly/2+0.5+dy
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt
    h['CRPIX3']=crpix
    h['CRVAL3']=crval
    h['CUNIT3']=('Angstrom','Units of coordinate increment and value    ')    
    h['CTYPE3']=('AWAV    ','Air wavelength (linear) ')
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='ICRS     '
    h['EQUINOX']=2000.00
    h['IFUCON']=(str(int(ns))+' ','NFibers')
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
        file=out_path+basenameC.replace('NAME',outf)
    out_fit=file
    hlist.writeto(out_fit,overwrite=True)
    tools.sycall('gzip -f '+out_fit)