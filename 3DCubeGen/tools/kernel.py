

def task_wrapper(args):
    return ifu_const(*args)
    
def ifu_const(spec_ifu,specE_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,npros,nproc,erroF):
#def ifu_const(spec_ifu0,specE_ifu0,x_ifu_V0,y_ifu_V0,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,it,npros,nproc,wcs,ra0,dec0):    
    val=int(nl/nproc)
    a1=val*npros
    if npros < nproc-1:
        a2=val*(npros+1)
    else:
        a2=nl
    spec_fint=[]
    specE_fint=[]
    if sigm_s > fibA*3.5*2:
        radiT=sigm_s/2.0
    else:
        radiT=fibA*3.5*2/2.0
    if len(x_ifu_V.shape) > 1:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            #Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            ##Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (fibA*3.5*2/2.0))[0]
            #spec_ifu,specE_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp,:],specE_ifu0[ntp,:],x_ifu_V0[ntp,:],y_ifu_V0[ntp,:]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spt_new=np.zeros(nw)
            if erroF:
                sptE_new=np.zeros(nw)
            Wgt=np.zeros(nw)
            Rspt=np.sqrt((y_ifu_V[:,0]-(yf+yi)/2.0)**2.0)
            #Rspt=np.sqrt((y_ifu_V[:,0]-yat)**2.0)
            ntpt=np.where(Rspt <= (radiT))[0]
            x_ifu_Vt=x_ifu_V[ntpt,:]
            y_ifu_Vt=y_ifu_V[ntpt,:]
            spec_ifut=spec_ifu[ntpt,:]
            if erroF:
                specE_ifut=specE_ifu[ntpt,:]
            for k in range(0, len(x_ifu_Vt[:,0])):
                Rsp=np.sqrt((x_ifu_Vt[k,:]-(xf+xi)/2.0)**2.0+(y_ifu_Vt[k,:]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_Vt[k,:]-xat)**2.0+(y_ifu_Vt[k,:]-yat)**2.0)
                ntp=np.where((Rsp <= (radiT)) & np.isfinite(spec_ifut[k,:]) & (spec_ifut[k,:] > 0))
                Wg=np.zeros(nw)
                if len(ntp[0]) > 0:   
                    Wg[ntp]=np.exp(-(Rsp[ntp]/sigm_s)**alph_s/2.0)
                    spt_new[ntp]=spec_ifut[k,ntp]*Wg[ntp]+spt_new[ntp]
                    if erroF:
                        sptE_new[ntp]=(specE_ifut[k,ntp]**2.0)*Wg[ntp]+sptE_new[ntp]
                Wgt=Wgt+Wg
            ntp=np.where(Wgt == 0)
            if len(ntp[0]) > 0:
                Wgt[ntp]=1
            spec_fint.extend([spt_new/Wgt])
            if erroF:
                specE_fint.extend([np.sqrt(sptE_new/Wgt)])
    else:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            #Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            ##Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (fibA*3.5*2/2.0))[0]
            #spec_ifu,specE_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp],specE_ifu0[ntp],x_ifu_V0[ntp],y_ifu_V0[ntp]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spax_new=0
            spaxE_new=0
            Wgt1=0
            Wgt2=0
            for k in range(0, len(x_ifu_V)):
                Rsp=np.sqrt((x_ifu_V[k]-(xf+xi)/2.0)**2.0+(y_ifu_V[k]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_V[k]-xat)**2.0+(y_ifu_V[k]-yat)**2.0)
                if Rsp <= (radiT):   
                    Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)
                    if np.isfinite(spec_ifu[k]): 
                        spax_new=spec_ifu[k]*Wg+spax_new
                        Wgt1=Wgt1+Wg
                    if np.isfinite(specE_ifu[k]):     
                        spaxE_new=(specE_ifu[k]**2.0)*Wg+spaxE_new
                        Wgt2=Wgt2+Wg
            if Wgt1 == 0:
                Wgt1=1
            if Wgt2 == 0:
                Wgt2=1    
            spec_fint.extend([spax_new/Wgt1])
            specE_fint.extend([np.sqrt(spaxE_new/Wgt2)])
    return ([spec_fint,specE_fint])