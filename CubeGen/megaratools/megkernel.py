import numpy as np


#spec_ifu,specE_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,npros,nproc,erroF
def kernel_int(spec_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yi,yf,xi,xf,nw,nl,erroF=False):
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
                ##fib=np.int(fib_idt[k])-1
                spt_new[ntp]=spec_ifu[k,ntp]*Wg[ntp]+spt_new[ntp]
                ##spt_err=(spec_ifu_e[k,:]*Wg)**2.0+spt_err**2.0
            Wgt=Wgt+Wg
        #if Wgt == 0:
        #    Wgt=1
        #print Wgt
        ntp=np.where(Wgt == 0)
        if len(ntp[0]) > 0:
            Wgt[ntp]=1
        ifu[:,j,i]=spt_new/Wgt
        ##if np.sum(np.sqrt(spt_err/Wgt**2.0)) == 0:
        ##    ifu_e[:,j,i]=1.0
        ##else:
        ##    ifu_e[:,j,i]=np.sqrt(spt_err/Wgt**2.0)
        #   # ifu_m[:,j,i]=1.0
        #int_spect=int_spect+spt_new/Wgt