import numpy as np

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
                        sptE_new[ntp]=(specE_ifut[k,ntp]**2.0)*Wg[ntp]**2.0+sptE_new[ntp]
                Wgt=Wgt+Wg
            ntp=np.where(Wgt == 0)
            if len(ntp[0]) > 0:
                Wgt[ntp]=1
            spec_fint.extend([spt_new/Wgt])
            if erroF:
                specE_fint.extend([np.sqrt(sptE_new)/Wgt])
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
                        spaxE_new=(specE_ifu[k]**2.0)*Wg**2.0+spaxE_new
                        Wgt2=Wgt2+Wg
            if Wgt1 == 0:
                Wgt1=1
            if Wgt2 == 0:
                Wgt2=1    
            spec_fint.extend([spax_new/Wgt1])
            specE_fint.extend([np.sqrt(spaxE_new)/Wgt2])
    return ([spec_fint,specE_fint])

def task_wrappermap(args):
    return map_const(*args)

def map_const(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,npros,nproc):    
#def ifu_const(spec_ifu0,specE_ifu0,specM_ifu0,specEM_ifu0,x_ifu_V0,y_ifu_V0,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,it,npros,nproc,wcs,ra0,dec0):    
    val=int(nl/nproc)
    a1=val*npros
    if npros < nproc-1:
        a2=val*(npros+1)
    else:
        a2=nl
    spec_fint=[]
    specE_fint=[]
    specM_fint=[]
    specEM_fint=[]
    if sigm_s > fibA*3.5*2:
        radiT=sigm_s/2.0
    else:
        radiT=fibA*3.5*2/2.0
    if len(x_ifu_V.shape) > 1:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            ##Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            #Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (radiT))[0]    
            #spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp,:],specE_ifu0[ntp,:],specM_ifu0[ntp,:],specEM_ifu0[ntp,:],x_ifu_V0[ntp,:],y_ifu_V0[ntp,:]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spt_new=np.zeros(nw)
            sptE_new=np.zeros(nw)
            sptM_new=np.zeros(nw)
            sptEM_new=np.zeros(nw)
            Wgt=np.zeros(nw)
            Rspt=np.sqrt((y_ifu_V[:,0]-(yf+yi)/2.0)**2.0)
            #Rspt=np.sqrt((y_ifu_V[:,0]-yat)**2.0)
            ntpt=np.where(Rspt <= (radiT))[0]
            x_ifu_Vt=x_ifu_V[ntpt,:]
            y_ifu_Vt=y_ifu_V[ntpt,:]
            spec_ifut=spec_ifu[ntpt,:]
            specE_ifut=specE_ifu[ntpt,:]
            specM_ifut=specM_ifu[ntpt,:]
            specEM_ifut=specEM_ifu[ntpt,:]
            for k in range(0, len(x_ifu_Vt[:,0])):
                Rsp=np.sqrt((x_ifu_Vt[k,:]-(xf+xi)/2.0)**2.0+(y_ifu_Vt[k,:]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_Vt[k,:]-xat)**2.0+(y_ifu_Vt[k,:]-yat)**2.0)
                ntp=np.where((Rsp <= (radiT)) & np.isfinite(spec_ifut[k,:]) & (spec_ifut[k,:] > 0))
                Wg=np.zeros(nw)
                if len(ntp[0]) > 0:   
                    Wg[ntp]=np.exp(-(Rsp[ntp]/sigm_s)**alph_s/2.0)
                    spt_new[ntp]=spec_ifut[k,ntp]*Wg[ntp]+spt_new[ntp]
                    sptE_new[ntp]=(specE_ifut[k,ntp]*Wg[ntp])**2.0+sptE_new[ntp]
                    sptM_new[ntp]=specM_ifut[k,ntp]*Wg[ntp]+sptM_new[ntp]
                    sptEM_new[ntp]=(specEM_ifut[k,ntp]*Wg[ntp])**2.0+sptEM_new[ntp]
                Wgt=Wgt+Wg
            ntp=np.where(Wgt == 0)
            if len(ntp[0]) > 0:
                Wgt[ntp]=1
            spec_fint.extend([spt_new/Wgt])
            specE_fint.extend([np.sqrt(sptE_new)/Wgt])
            specM_fint.extend([sptM_new/Wgt])
            specEM_fint.extend([np.sqrt(sptEM_new)/Wgt])
    else:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            ##Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            #Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (radiT))[0]
            #spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp],specE_ifu0[ntp],specM_ifu0[ntp],specEM_ifu0[ntp],x_ifu_V0[ntp],y_ifu_V0[ntp]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spax_new=0
            spaxE_new=0
            spaxM_new=0
            spaxEM_new=0
            Wgt1=0
            Wgt2=0
            Wgt3=0
            Wgt4=0
            for k in range(0, len(x_ifu_V)):
                Rsp=np.sqrt((x_ifu_V[k]-(xf+xi)/2.0)**2.0+(y_ifu_V[k]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_V[k]-xat)**2.0+(y_ifu_V[k]-yat)**2.0)
                if Rsp <= (radiT):   
                    Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)
                    if np.isfinite(spec_ifu[k]):
                        spax_new=spec_ifu[k]*Wg+spax_new
                        Wgt1=Wgt1+Wg
                    if np.isfinite(specE_ifu[k]):
                        spaxE_new=(specE_ifu[k]*Wg)**2.0+spaxE_new
                        Wgt2=Wgt2+Wg
                    if np.isfinite(specM_ifu[k]):    
                        spaxM_new=specM_ifu[k]*Wg+spaxM_new
                        Wgt3=Wgt3+Wg
                    if np.isfinite(specEM_ifu[k]):    
                        spaxEM_new=(specEM_ifu[k]*Wg)**2.0+spaxEM_new
                        Wgt4=Wgt4+Wg
            if Wgt1 == 0:
                Wgt1=1
            if Wgt2 == 0:
                Wgt2=1
            if Wgt3 == 0:
                Wgt3=1
            if Wgt4 == 0:
                Wgt4=1    
            spec_fint.extend([spax_new/Wgt1])
            specE_fint.extend([np.sqrt(spaxE_new)/Wgt2])
            specM_fint.extend([spaxM_new/Wgt3])
            specEM_fint.extend([np.sqrt(spaxEM_new)/Wgt4])
    return ([spec_fint,specE_fint,specM_fint,specEM_fint])

def task_wrappermatrix(args):
    return matrix_const(*args)


def matrix_const(spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,npros,nproc):
#def ifu_const(spec_ifu0,specE_ifu0,specM_ifu0,specEM_ifu0,x_ifu_V0,y_ifu_V0,fibA,pix_s,sigm_s,alph_s,yo,xi,xf,nw,nl,it,npros,nproc,wcs,ra0,dec0):    
    val=int(nl/nproc)
    a1=val*npros
    if npros < nproc-1:
        a2=val*(npros+1)
    else:
        a2=nl
    spec_fint=[]
    specE_fint=[]
    specM_fint=[]
    specEM_fint=[]
    if sigm_s > fibA*3.5*2:
        radiT=sigm_s/2.0
    else:
        radiT=fibA*3.5*2/2.0
    if len(x_ifu_V.shape) > 1:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            ##Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            #Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (radiT))[0]    
            #spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp,:],specE_ifu0[ntp,:],specM_ifu0[ntp,:],specEM_ifu0[ntp,:],x_ifu_V0[ntp,:],y_ifu_V0[ntp,:]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spt_new=np.zeros(nw)
            sptE_new=np.zeros(nw)
            sptM_new=np.zeros(nw)
            sptEM_new=np.zeros(nw)
            Wgt=np.zeros(nw)
            Rspt=np.sqrt((y_ifu_V[:,0]-(yf+yi)/2.0)**2.0)
            #Rspt=np.sqrt((y_ifu_V[:,0]-yat)**2.0)
            ntpt=np.where(Rspt <= (radiT))[0]
            x_ifu_Vt=x_ifu_V[ntpt,:]
            y_ifu_Vt=y_ifu_V[ntpt,:]
            spec_ifut=spec_ifu[ntpt,:]
            specE_ifut=specE_ifu[ntpt,:]
            specM_ifut=specM_ifu[ntpt,:]
            specEM_ifut=specEM_ifu[ntpt,:]
            for k in range(0, len(x_ifu_Vt[:,0])):
                Rsp=np.sqrt((x_ifu_Vt[k,:]-(xf+xi)/2.0)**2.0+(y_ifu_Vt[k,:]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_Vt[k,:]-xat)**2.0+(y_ifu_Vt[k,:]-yat)**2.0)
                ntp=np.where((Rsp <= (radiT)) & np.isfinite(spec_ifut[k,:]) & (spec_ifut[k,:] > 0))
                Wg=np.zeros(nw)
                if len(ntp[0]) > 0:   
                    Wg[ntp]=np.exp(-(Rsp[ntp]/sigm_s)**alph_s/2.0)
                    spt_new[ntp]=spec_ifut[k,ntp]*Wg[ntp]+spt_new[ntp]
                    sptE_new[ntp]=(specE_ifut[k,ntp]*Wg[ntp])**2.0+sptE_new[ntp]
                    sptM_new[ntp]=specM_ifut[k,ntp]*Wg[ntp]+sptM_new[ntp]
                    sptEM_new[ntp]=(specEM_ifut[k,ntp]*Wg[ntp])**2.0+sptEM_new[ntp]
                Wgt=Wgt+Wg
            ntp=np.where(Wgt == 0)
            if len(ntp[0]) > 0:
                Wgt[ntp]=1
            spec_fint.extend([spt_new/Wgt])
            specE_fint.extend([np.sqrt(sptE_new)/Wgt])
            specM_fint.extend([sptM_new/Wgt])
            specEM_fint.extend([np.sqrt(sptEM_new)/Wgt])
    else:
        for j in range(a1, a2):
            #skycor = pixel_to_skycoord(it,j,wcs)
            #xat=-(skycor.ra.value*3600.0-ra0)
            #yat=skycor.dec.value*3600.0-dec0
            ##Rsp=np.sqrt((x_ifu_V0-(xf+xi)/2.0)**2.0)
            #Rsp=np.sqrt((x_ifu_V0-xat)**2.0)
            #ntp=np.where(Rsp <= (radiT))[0]
            #spec_ifu,specE_ifu,specM_ifu,specEM_ifu,x_ifu_V,y_ifu_V=spec_ifu0[ntp],specE_ifu0[ntp],specM_ifu0[ntp],specEM_ifu0[ntp],x_ifu_V0[ntp],y_ifu_V0[ntp]
            yi=yo+pix_s*j
            yf=yo+pix_s*(j+1)
            spax_new=0
            spaxE_new=0
            spaxM_new=0
            spaxEM_new=0
            Wgt1=0
            Wgt2=0
            Wgt3=0
            Wgt4=0
            for k in range(0, len(x_ifu_V)):
                Rsp=np.sqrt((x_ifu_V[k]-(xf+xi)/2.0)**2.0+(y_ifu_V[k]-(yf+yi)/2.0)**2.0)
                #Rsp=np.sqrt((x_ifu_V[k]-xat)**2.0+(y_ifu_V[k]-yat)**2.0)
                if Rsp <= (radiT):   
                    Wg=np.exp(-(Rsp/sigm_s)**alph_s/2.0)
                    if np.isfinite(spec_ifu[k]):
                        spax_new=spec_ifu[k]*Wg+spax_new
                        Wgt1=Wgt1+Wg
                    if np.isfinite(specE_ifu[k]):
                        spaxE_new=(specE_ifu[k]*Wg)**2.0+spaxE_new
                        Wgt2=Wgt2+Wg
                    if np.isfinite(specM_ifu[k]):    
                        spaxM_new=specM_ifu[k]*Wg+spaxM_new
                        Wgt3=Wgt3+Wg
                    if np.isfinite(specEM_ifu[k]):    
                        spaxEM_new=(specEM_ifu[k]*Wg)**2.0+spaxEM_new
                        Wgt4=Wgt4+Wg
            if Wgt1 == 0:
                Wgt1=1
            if Wgt2 == 0:
                Wgt2=1
            if Wgt3 == 0:
                Wgt3=1
            if Wgt4 == 0:
                Wgt4=1    
            spec_fint.extend([spax_new/Wgt1])
            specE_fint.extend([np.sqrt(spaxE_new)/Wgt2])
            specM_fint.extend([spaxM_new/Wgt3])
            specEM_fint.extend([np.sqrt(spaxEM_new)/Wgt4])
    return ([spec_fint,specE_fint,specM_fint,specEM_fint])


    def task_wrappercov1(args):
        return covdotone(*args)

    def covdotone(St,Wg,nx,ns,a1,a2):
        dy=a2-a1
        outt=np.zeros([nx,dy,ns])
        ct=0
        for i in range(a1, a2):
            t=np.dot(St**2, Wg[:,ct,:].transpose(1,0))
            outt[:,ct,:]=t.transpose(1,0)
            ct=ct+1
        return ([outt])

    def task_wrappercov2(args):
        return covdottwo(*args)

    def covdottwo(St,Wg,out,Dq,nx,i,a1,a2):
        dy=a2-a1
        outt2=np.zeros([nx,nx,dy])
        distF=np.zeros([nx,nx,dy])
        ct=0
        for j in range(a1, a2):
            dt=np.sqrt(squareform(Dq)**2+(i-j)**2)
            t=np.dot(Wg,out[:,ct,:].transpose(1,0))
            outt2[:,:,ct]=t.transpose(1,0)
            distF[:,:,ct]=dt
            ct=ct+1
        return ([outt2,distF])