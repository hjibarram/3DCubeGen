import numpy as np


#spec_ifu,specE_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,npros,nproc,erroF
def kernel_int(ifu,ifuE,spec_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yi,yf,xi,xf,nw,nl,i,erroF=False):
    if sigm_s > fibA*3.5*2:
        radiT=sigm_s/2.0
    else:
        radiT=fibA*3.5/2.0
    for j in range(0, nl):
        yi=yf
        yf=yf+pix_s
        spt_new=np.zeros(nw)
        if erroF:
            sptE_new=np.zeros(nw)
        Wgt=np.zeros(nw)
        for k in range(0, len(x_ifu_V[:,0])):
            Rsp=np.sqrt((x_ifu_V[k,:]-(xf+xi)/2.0)**2.0+(y_ifu_V[k,:]-(yf+yi)/2.0)**2.0)
            ntp=np.where((Rsp <= (radiT)) & np.isfinite(spec_ifu[k,:]) & (spec_ifu[k,:] > 0))
                
            Wg=np.zeros(nw)
            if len(ntp[0]) > 0:   
                Wg[ntp]=np.exp(-(Rsp[ntp]/sigm_s)**alph_s/2.0)
                spt_new[ntp]=spec_ifu[k,ntp]*Wg[ntp]+spt_new[ntp]
                if erroF:
                    sptE_new[ntp]=(specE_ifut[k,ntp]**2.0)*Wg[ntp]**2.0+sptE_new[ntp]
            Wgt=Wgt+Wg
        ntp=np.where(Wgt == 0)
        if len(ntp[0]) > 0:
            Wgt[ntp]=1
        ifu[:,j,i]=spt_new/Wgt
        if erroF:
            ifuE[:,j,i]=np.sqrt(sptE_new)/Wgt
    return ifu,ifuE


def kernel_intmap(ifu,ifuE,spec_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yi,yf,xi,xf,nl,i,erroF=False):
    if sigm_s > fibA*3.5*2:
        radiT=sigm_s/2.0
    else:
        radiT=fibA*3.5/2.0
    for j in range(0, nl):
        yi=yf
        yf=yf+pix_s
        spax_new=0
        Wgt1=0
        if erroF:
            Wgt2=0
            spaxE_new=0    
        for k in range(0, len(x_ifu_V)):
            Rsp=np.sqrt((x_ifu_V[k]-(xf+xi)/2.0)**2.0+(y_ifu_V[k]-(yf+yi)/2.0)**2.0)
            if Rsp <= (radiT):   
                Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)
                if np.isfinite(spec_ifu[k]):
                    spax_new=spec_ifu[k]*Wg+spax_new
                    Wgt1=Wgt1+Wg
                if erroF:
                    if np.isfinite(specE_ifu[k]):
                        spaxE_new=(specE_ifu[k]*Wg)**2.0+spaxE_new
                        Wgt2=Wgt2+Wg
        if Wgt1 == 0:
            Wgt1=1 
        ifu[j,i]=spax_new/Wgt1
        if erroF:
            if Wgt2 == 0:
                Wgt2=1       
            ifuE[j,i]=np.sqrt(spaxE_new)/Wgt2
    return ifu,ifuE