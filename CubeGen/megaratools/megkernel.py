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
        for k in range(0, len(x_ifu)):
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