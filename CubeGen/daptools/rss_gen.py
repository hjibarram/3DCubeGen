from tqdm.notebook import tqdm
from tqdm import tqdm as tqdmT
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import CubeGen.tools.tools as tools
from astropy.wcs.utils import pixel_to_skycoord

def rssp_extract(name,path='./',path_out='./',basename_in='lvmCube-NAME.fits.gz',basename_out='lvmRSS-NAMElab.fits',flu16=True,nsplit=0,spt=[0,0],lvm=True,fluxu=1e-16,notebook=True):
    file=basename_in.replace('NAME',name)
    cube_file=path+'/'+file
    [flux, hdr]=fits.getdata(cube_file, 0, header=True)
    if flu16:
        flux=flux*fluxu
    try:
        if lvm:
            hdrO=fits.getheader(cube_file, 3)
            orhd=True
        else:
            orhd=False
    except:
        orhd=False

    crpix=hdr["CRPIX3"]
    try:
        cdelt=hdr["CD3_3"]
    except:
        cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    nz,nx0,ny0=flux.shape
    wave=crval+cdelt*(np.arange(nz)+1-crpix)
    wcs=WCS(hdr)
    wcs=wcs.celestial

    if orhd:
        try:
            print('Loaded original header')
            nxo0=hdrO['NX0']
            nyo0=hdrO['NY0']
            nxo1=hdrO['NX1']
            nyo1=hdrO['NY1']
            nxo2=hdrO['NX2']
            nyo2=hdrO['NY2']
        except:
            nxo0=nx0
            nyo0=ny0
            nxo1=0
            nyo1=0
            nxo2=nx0
            nyo2=ny0
    else:
        nxo0=nx0
        nyo0=ny0
        nxo1=0
        nyo1=0
        nxo2=nx0
        nyo2=ny0
    
    if nsplit > 1:
        nx_list=[]
        ny_list=[]
        dxt=int(nx0/nsplit)
        dyt=int(ny0/nsplit)
        for kt in range(0, nsplit):
            nx_list.extend([dxt*kt])
            ny_list.extend([dyt*kt])
        nx_list.extend([nx0])
        ny_list.extend([ny0])
        nx=nx_list[spt[0]+1]
        ny=ny_list[spt[1]+1]
        nx1=nx_list[spt[0]]
        ny1=ny_list[spt[1]]
        labf='_p'+str(spt[0])+str(spt[1])
        cubeT=np.copy(flux[:,nx1:nx,ny1:ny])
        del flux
        flux=np.copy(cubeT)
        del cubeT
    else:
        nx=nx0
        ny=ny0
        nx1=0
        ny1=0
        labf=''
    [fluxE, hdre]=fits.getdata(cube_file, 1, header=True)
    if nsplit > 1:
        cubeT=np.copy(fluxE[:,nx1:nx,ny1:ny])
        del fluxE
        fluxE=np.copy(cubeT)
        del cubeT
    if flu16:
        fluxE=fluxE*fluxu
    if notebook:
        pbar=tqdm(total=(nx-nx1)*(ny-ny1))
    else:
        pbar=tqdmT(total=(nx-nx1)*(ny-ny1))
    ct=0
    for i in range(nx1, nx):
        for j in range(ny1, ny):
            spec=flux[:,i-nx1,j-ny1]
            if np.nansum(spec) != 0:
                ct=ct+1
            pbar.update(1)
    pbar.close()
    ns=ct
    rss=np.zeros([ns,nz])
    rssI=np.zeros([ns,nz])
    rssM=np.zeros([ns,nz])
    rssW=np.zeros([1,nz])
    rssL=np.ones([ns,nz])
    rssS1=np.ones([ns,nz])
    rssS2=np.ones([ns,nz])
    rssS3=np.ones([ns,nz])
    rssS4=np.ones([ns,nz])
    fiberid=np.zeros(ns)
    specid=np.ones(ns)
    blockid=[]
    finblock=np.zeros(ns)
    targettype=[]
    ifulabel=[]
    finifu=np.zeros(ns)
    telescope=[]
    xpmm=np.zeros(ns)
    ypmm=np.zeros(ns)
    ringnum=np.ones(ns)
    orig_ifulabel=[]
    orig_slitlabel=[]
    finsector=np.zeros(ns)
    fmap=[]
    ypix_b=np.zeros(ns)
    ypix_r=np.zeros(ns)
    ypix_z=np.zeros(ns)
    fibstatus=np.zeros(ns)
    ra=np.zeros(ns)
    dec=np.zeros(ns)
    rssW[0,:]=wave
    ct=0
    if notebook:
        pbar=tqdm(total=ns)
    else:
        pbar=tqdmT(total=ns)
    for i in range(nx1, nx):
        for j in range(ny1, ny):
            spec=flux[:,i-nx1,j-ny1]
            specE=fluxE[:,i-nx1,j-ny1]
            if np.nansum(spec) != 0:
                rss[ct,:]=spec
                rssI[ct,:]=1/specE**2.0
                fiberid[ct]=ct+1
                blockid.extend(['B1'])
                finblock[ct]=ct+1
                targettype.extend(['science'])
                ifulabel.extend(['Sci1'])
                finifu[ct]=ct+1
                telescope.extend(['Sci'])
                xpmm[ct]=j+nyo1
                ypmm[ct]=i+nxo1
                orig_ifulabel.extend(['S1-'+str(ct+1)])
                orig_slitlabel.extend(['S1B1-'+str(ct+1)])
                finsector[ct]=ct+1
                fmap.extend(['S1-'+str(ct+1)+':'+'S1B1-'+str(ct+1)])
                sky1=pixel_to_skycoord(j,i,wcs)
                ra[ct]=sky1.ra.value
                dec[ct]=sky1.dec.value
                ct=ct+1
                pbar.update(1)
    pbar.close()
    blockid=np.array(blockid)
    targettype=np.array(targettype)
    ifulabel=np.array(ifulabel)
    telescope=np.array(telescope)
    orig_ifulabel=np.array(orig_ifulabel)
    orig_slitlabel=np.array(orig_slitlabel)
    fmap=np.array(fmap)
    
    tools.sycall('mkdir -p '+path_out)

    h1=fits.PrimaryHDU()
    h2=fits.ImageHDU(rss)
    h3=fits.ImageHDU(rssI)
    h4=fits.ImageHDU(rssM)
    h5=fits.ImageHDU(rssW)
    h6=fits.ImageHDU(rssL)
    h7=fits.ImageHDU(rssS1)
    h8=fits.ImageHDU(rssS2)
    h9=fits.ImageHDU(rssS3)
    h10=fits.ImageHDU(rssS4)
    h_k=h1.header
    keys=list(hdr.keys())
    for i in range(0, len(keys)):
        h_k[keys[i]]=hdr[keys[i]]
        h_k.comments[keys[i]]=hdr.comments[keys[i]]
    h_k['EXTNAME']='ORIGINAL'
    h_k['POSCIRA']=hdr['CRVAL1']                                        
    h_k['POSCIDE']=hdr['CRVAL2']   
    h_k['EXPOSURE']=000000
    h_k.update()
    h_t=h2.header
    h_t['EXTNAME']='FLUX'
    h_t['CDELT1']=cdelt
    h_t['CRPIX1']=crpix
    h_t['CRVAL1']=crval
    h_t['CUNIT1']='Angstrom'  
    h_t['CTYPE1']='WAVE    '
    h_t['CDELT2']=1
    h_t['CRPIX2']=1
    h_t['CRVAL2']=1
    h_t['CUNIT2']=''  
    h_t['CTYPE2']='FIBERID '
    h_t['BUNIT']='erg/s/cm^2'
    h_t.update()
    h_r=h3.header
    h_r['EXTNAME']='IVAR'
    h_r['CDELT1']=cdelt
    h_r['CRPIX1']=crpix
    h_r['CRVAL1']=crval
    h_r['CUNIT1']='Angstrom'  
    h_r['CTYPE1']='WAVE    '
    h_r['CDELT2']=1
    h_r['CRPIX2']=1
    h_r['CRVAL2']=1
    h_r['CUNIT2']=''  
    h_r['CTYPE2']='FIBERID '
    h_r['BUNIT']='erg/s/cm^2'
    h_r.update()    
    h_m=h4.header
    h_m['EXTNAME']='MASK'
    h_m['CDELT1']=cdelt
    h_m['CRPIX1']=crpix
    h_m['CRVAL1']=crval
    h_m['CUNIT1']='Angstrom'  
    h_m['CTYPE1']='WAVE    '
    h_m['CDELT2']=1
    h_m['CRPIX2']=1
    h_m['CRVAL2']=1
    h_m['CUNIT2']=''  
    h_m['CTYPE2']='FIBERID '
    h_m['BUNIT']='erg/s/cm^2'
    h_m.update()    
    h_w=h5.header
    h_w['EXTNAME']='WAVE'
    h_w['CDELT1']=cdelt
    h_w['CRPIX1']=crpix
    h_w['CRVAL1']=crval
    h_w['CUNIT1']='Angstrom'  
    h_w['CTYPE1']='WAVE    '
    h_w['CDELT2']=1
    h_w['CRPIX2']=1
    h_w['CRVAL2']=1
    h_w['CUNIT2']=''  
    h_w['CTYPE2']='FIBERID '
    h_w.update()    
    h_l=h6.header
    h_l['EXTNAME']='LSF'
    h_l['CDELT1']=cdelt
    h_l['CRPIX1']=crpix
    h_l['CRVAL1']=crval
    h_l['CUNIT1']='Angstrom'  
    h_l['CTYPE1']='WAVE    '
    h_l['CDELT2']=1
    h_l['CRPIX2']=1
    h_l['CRVAL2']=1
    h_l['CUNIT2']=''  
    h_l['CTYPE2']='FIBERID '
    #COMMENT LSF (FWHM) solution in angstroms
    h_l.update()
    h_s=h7.header
    h_s['EXTNAME']='SKY_EAST'
    h_s['CDELT1']=cdelt
    h_s['CRPIX1']=crpix
    h_s['CRVAL1']=crval
    h_s['CUNIT1']='Angstrom'  
    h_s['CTYPE1']='WAVE    '
    h_s['CDELT2']=1
    h_s['CRPIX2']=1
    h_s['CRVAL2']=1
    h_s['CUNIT2']=''  
    h_s['CTYPE2']='FIBERID '
    #COMMENT sky east in flux-calibrated units (10(-17) erg/s/cm2/Ang/fiber) 
    h_s.update()
    h_s1=h8.header
    h_s1['EXTNAME']='SKY_EAST_IVAR'
    h_s1['CDELT1']=cdelt
    h_s1['CRPIX1']=crpix
    h_s1['CRVAL1']=crval
    h_s1['CUNIT1']='Angstrom'  
    h_s1['CTYPE1']='WAVE    '
    h_s1['CDELT2']=1
    h_s1['CRPIX2']=1
    h_s1['CRVAL2']=1
    h_s1['CUNIT2']=''  
    h_s1['CTYPE2']='FIBERID '
    #COMMENT sky east in flux-calibrated units (10(-17) erg/s/cm2/Ang/fiber) 
    h_s1.update()        
    h_s2=h9.header
    h_s2['EXTNAME']='SKY_WEST'
    h_s2['CDELT1']=cdelt
    h_s2['CRPIX1']=crpix
    h_s2['CRVAL1']=crval
    h_s2['CUNIT1']='Angstrom'  
    h_s2['CTYPE1']='WAVE    '
    h_s2['CDELT2']=1
    h_s2['CRPIX2']=1
    h_s2['CRVAL2']=1
    h_s2['CUNIT2']=''  
    h_s2['CTYPE2']='FIBERID '
    #COMMENT sky east in flux-calibrated units (10(-17) erg/s/cm2/Ang/fiber) 
    h_s2.update()
    h_s3=h10.header
    h_s3['EXTNAME']='SKY_WEST_IVAR'
    h_s3['CDELT1']=cdelt
    h_s3['CRPIX1']=crpix
    h_s3['CRVAL1']=crval
    h_s3['CUNIT1']='Angstrom'  
    h_s3['CTYPE1']='WAVE    '
    h_s3['CDELT2']=1
    h_s3['CRPIX2']=1
    h_s3['CRVAL2']=1
    h_s3['CUNIT2']=''  
    h_s3['CTYPE2']='FIBERID '
    #COMMENT sky east in flux-calibrated units (10(-17) erg/s/cm2/Ang/fiber) 
    h_s3.update() 
    col1 = fits.Column(name='fiberid', format='K', array=fiberid)
    col2 = fits.Column(name='spectrographid', format='K', array=specid)
    col3 = fits.Column(name='blockid', format='3A', array=blockid)
    col4 = fits.Column(name='finblock', format='K', array=finblock)
    col5 = fits.Column(name='targettype', format='8A', array=targettype)
    col6 = fits.Column(name='ifulabel', format='5A', array=ifulabel)
    col7 = fits.Column(name='finifu', format='K', array=finifu)
    col8 = fits.Column(name='telescope', format='4A', array=telescope)
    col9 = fits.Column(name='xpmm', format='D', array=xpmm)
    col10 = fits.Column(name='ypmm', format='D', array=ypmm)
    col11 = fits.Column(name='ringnum', format='D', array=ringnum)
    col12 = fits.Column(name='orig_ifulabel', format='6A', array=orig_ifulabel)
    col13 = fits.Column(name='orig_slitlabel', format='8A', array=orig_slitlabel)
    col14 = fits.Column(name='finsector', format='K', array=finsector)
    col15 = fits.Column(name='fmap', format='17A', array=fmap)
    col16 = fits.Column(name='ypix_b', format='K', array=ypix_b)
    col17 = fits.Column(name='ypix_r', format='K', array=ypix_r)
    col18 = fits.Column(name='ypix_z', format='K', array=ypix_z)
    col19 = fits.Column(name='fibstatus', format='K', array=fibstatus)
    col20 = fits.Column(name='ra', format='D', array=ra)
    col21 = fits.Column(name='dec', format='D', array=dec)
    coldefs = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                            col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21])
    h11 = fits.BinTableHDU.from_columns(coldefs)
    h_y=h11.header
    h_y['EXTNAME']='SLITMAP'
    h_y.update()
    hlist=fits.HDUList([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11])
    hlist.update_extend()
    hlist.writeto(path_out+'/'+basename_out.replace(labf,'').replace('NAME',name).replace('lab',labf), overwrite=True)
    tools.sycall('gzip -f '+path_out+'/'+basename_out.replace(labf,'').replace('NAME',name).replace('lab',labf))