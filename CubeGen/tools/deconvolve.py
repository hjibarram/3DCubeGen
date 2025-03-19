import matplotlib.pyplot as plt
from sdeconv.data import celegans
from sdeconv.psfs import SPSFGaussian
from sdeconv.deconv import SRichardsonLucy
from astropy.io import fits
import torch
import numpy as np
from skimage import color, data, restoration
from scipy.signal import convolve2d as conv2
from tqdm.notebook import tqdm
import CubeGen.tools.tools as tools


def psfG(psf_x=1.33,psf_y=1.33,nx=35,ny=35):
    #psf_s=1.33*1.5
    # Generate a 2D PSF
    psf_generator = SPSFGaussian((psf_x, psf_y), (nx, ny))
    psf = psf_generator()
    #psf.numpy()
    return psf

def deconvolve_2dmap(image,psf_x=1.33,psf_y=1.33,nxpsf=35,nypsf=35,niter=10):
    psf=psfG(psf_x=psf_x,psf_y=psf_y,nx=nxpsf,ny=nypsf)
    imageT=np.copy(image)
    imageT=imageT/np.nanmax(imageT)
    deconvolved_RL = restoration.richardson_lucy(imageT, psf, num_iter=niter)
    fac=np.nanmean(image)/np.nanmean(deconvolved_RL)
    deconvolved_RL=deconvolved_RL*fac
    return deconvolved_RL

def deconvolve_3dcube(name='name',path='',psf_x=1.33,psf_y=1.33,nxpsf=35,nypsf=35,niter=10,pbars=True,basenameC='lvmCube-NAME.fits'):
    file=path+basenameC.replace('NAME',name)
    try:
        [cube,hdr1]=fits.getdata(file, 0, header=True)
    except:
        file=file+'.gz'
        [cube,hdr1]=fits.getdata(file, 0, header=True)
    try:
        [cubeE,hdr2]=fits.getdata(file, 1, header=True)
        error=True
    except:
        error=False
    nz,nx,ny=cube.shape
    cubeT=np.copy(cube)
    if pbars:    
        pbar=tqdm(total=nz) 
    for i in range(0, nz):
        image=cube[i,:,:]
        imageT=deconvolve_2dmap(image,psf_x=psf_x,psf_y=psf_y,nxpsf=nxpsf,nypsf=nypsf,niter=niter)
        cubeT[i,:,:]=imageT
        if pbars:
            pbar.update(1)     
    if pbars:
        pbar.close()

    outfile1=path+basenameC.replace('NAME',name+'_decv') 
    h1=fits.PrimaryHDU(cubeT,header=hdr1)
    if error:
        h2=fits.ImageHDU(cubeE,header=hdr2)
        hlist=fits.HDUList([h1,h2])
    else:
        hlist=fits.HDUList([h1])
    hlist.update_extend()
    hlist.writeto(outfile1,overwrite=True)
    tools.sycall('gzip -f '+outfile1)


def deconvolve_2dfile(name='name',path='',psf_x=1.33,psf_y=1.33,nxpsf=35,nypsf=35,niter=10,pbars=True,basenameC='lvmMap-NAME.fits'):
    file=path+basenameC.replace('NAME',name)
    try:
        [mapt,hdr1]=fits.getdata(file, 0, header=True)
    except:
        file=file+'.gz'
        [mapt,hdr1]=fits.getdata(file, 0, header=True)
    try:
        [maptM,hdr2]=fits.getdata(file, 1, header=True)
        magt=True
    except:
        magt=False
    if magt:
        try:
            [maptE,hdr3]=fits.getdata(file, 2, header=True)
            error=True
        except:
            error=False
    else:
        try:
            [maptE,hdr2]=fits.getdata(file, 1, header=True)
            error=True
        except:
            error=False
    maptT=deconvolve_2dmap(mapt,psf_x=psf_x,psf_y=psf_y,nxpsf=nxpsf,nypsf=nypsf,niter=niter)
    if magt:
        maptMT=deconvolve_2dmap(maptM,psf_x=psf_x,psf_y=psf_y,nxpsf=nxpsf,nypsf=nypsf,niter=niter)
  
    outfile1=path+basenameC.replace('NAME',name+'_decv') 
    h1=fits.PrimaryHDU(maptT,header=hdr1)
    if mag:
        h2=fits.ImageHDU(maptMT,header=hdr2)
        if error:
            h3=fits.ImageHDU(maptE,header=hdr3)
            hlist=fits.HDUList([h1,h2,h3])
        else:
            hlist=fits.HDUList([h1,h2])
    else:
        if error:
            h2=fits.ImageHDU(maptE,header=hdr2)
            hlist=fits.HDUList([h1,h2])
        else:
            hlist=fits.HDUList([h1])
    hlist.update_extend()
    hlist.writeto(outfile1,overwrite=True)
    tools.sycall('gzip -f '+outfile1)