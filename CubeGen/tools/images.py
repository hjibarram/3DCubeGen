from astropy.io import fits
import numpy as np
from PIL import Image as im2

def get_jpgNR(namel=['name1','name2','name3'],dir1='',dir0='',valTs=[[7,15],[7,13.5],[7,15.5]],name='image_rgb',wgt=[1.0,1.0,1.0]):
    [pdl_Ba,hdrt]=fits.getdata(dir0+'/'+namel[2], 'MAG', header=True)
    [pdl_Ga,hdrt]=fits.getdata(dir0+'/'+namel[1], 'MAG', header=True)
    [pdl_Ra,hdrt]=fits.getdata(dir0+'/'+namel[0], 'MAG', header=True)
    min1=valTs[2][1]
    min2=valTs[1][1]
    min3=valTs[0][1]
    max1=valTs[2][0]
    max2=valTs[1][0]
    max3=valTs[0][0]
    nx,ny=pdl_Ba.shape
    pdl_Ba=(np.flipud(pdl_Ba)-min1)/(max1-min1)*256
    pdl_Ga=(np.flipud(pdl_Ga)-min2)/(max2-min2)*256
    pdl_Ra=(np.flipud(pdl_Ra)-min3)/(max3-min3)*256

    pdl_B=wgt[2]*pdl_Ba
    pdl_G=wgt[1]*pdl_Ga
    pdl_R=wgt[0]*pdl_Ra
    #print(wgt[2])
    if len(namel) > 3:
        if len(namel) % 3 == 0:
            wt2=wgt[2]
            wt1=wgt[1]
            wt0=wgt[0]
            nt=int(len(namel)/3)-1
            for i in range(1, nt+1):
                [pdl_Bt,hdrt]=fits.getdata(dir0+'/'+namel[i*3+2], 'MAG', header=True)
                [pdl_Gt,hdrt]=fits.getdata(dir0+'/'+namel[i*3+1], 'MAG', header=True)
                [pdl_Rt,hdrt]=fits.getdata(dir0+'/'+namel[i*3+0], 'MAG', header=True)
                min1=valTs[i*3+2][1]
                min2=valTs[i*3+1][1]
                min3=valTs[i*3+0][1]
                max1=valTs[i*3+2][0]
                max2=valTs[i*3+1][0]
                max3=valTs[i*3+0][0]
                pdl_Bt=(np.flipud(pdl_Bt)-min1)/(max1-min1)*256
                pdl_Gt=(np.flipud(pdl_Gt)-min2)/(max2-min2)*256
                pdl_Rt=(np.flipud(pdl_Rt)-min3)/(max3-min3)*256
                pdl_B=pdl_B+wgt[i*3+2]*pdl_Bt
                pdl_G=pdl_G+wgt[i*3+1]*pdl_Gt
                pdl_R=pdl_R+wgt[i*3+0]*pdl_Rt
                wt2=wt2+wgt[i*3+2]
                wt1=wt1+wgt[i*3+1]
                wt0=wt0+wgt[i*3+0]
            pdl_B=pdl_B/wt2
            pdl_G=pdl_G/wt1
            pdl_R=pdl_R/wt0      
    #pdl_B=(pdl_B-np.min(pdl_B))/(np.max(pdl_B)-np.min(pdl_B))*256
    #pdl_G=(pdl_G-np.min(pdl_G))/(np.max(pdl_G)-np.min(pdl_G))*256
    #pdl_R=(pdl_R-np.min(pdl_R))/(np.max(pdl_R)-np.min(pdl_R))*256
    pdl_B[np.where(pdl_B < 0)]=0
    pdl_G[np.where(pdl_G < 0)]=0
    pdl_R[np.where(pdl_R < 0)]=0
    pdl_B[np.where(pdl_B > 255)]=255
    pdl_G[np.where(pdl_G > 255)]=255
    pdl_R[np.where(pdl_R > 255)]=255
    pdl_00=np.zeros([nx,ny,3],dtype="uint8")
    pdl_00[:,:,0]=pdl_R
    pdl_00[:,:,1]=pdl_G
    pdl_00[:,:,2]=pdl_B
    im1 = im2.fromarray(pdl_00)
    im1.save(dir1+'/'+name+'.jpeg',quality=100)