import numpy as np
from scipy.interpolate import interp1d

def extintion_c(wave,dir_tem='data',basename='extintion_curve'):
    f=open(dir_tem+'/'+basename+'.txt','r')
    x1=[]
    y1=[]
    for line in f:
        data=line.replace('\n','').split(' ')
        data=list(filter(None,data))
        x1.extend([float(data[0])])
        y1.extend([float(data[1])])
    y1=np.array(y1)
    x1=np.array(x1)
    Kl=interp1d(x1,y1,bounds_error=False,fill_value=0.)(wave)
    return Kl

def id_str(id,n_z=2):
    id=int(float(id))
    if n_z < 2 or n_z > 8:
        n_z=2
    if n_z == 2:
        if id < 10:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 3:
        if id < 10:
            idt='00'+str(id)
        elif id < 100:
            idt='0'+str(id)
        else:
            idt=str(id)
    return idt

def megarafiber_pos(hdr):
    nfib=hdr['NFIBERS']
    psc=hdr['PSCALE']
    x_pos=np.zeros(nfib)
    y_pos=np.zeros(nfib)
    fib_a=np.zeros(nfib)
    fib_b=np.zeros(nfib)
    fib_id=np.zeros(nfib)
    for i in range(0, nfib):
        #print 'FIB'+str(i+1)+'_X'
        x_pos[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_X']
        y_pos[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_Y']
        fib_a[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_D']
        fib_b[i]=hdr['FIB'+id_str(i+1,n_z=3)+'_B']
        fib_id[i]=i+1
    nt=np.where((np.abs(x_pos) < 10) & (np.abs(y_pos) < 10))
    x_posf=x_pos[nt]
    y_posf=y_pos[nt]
    fib_idt=fib_id[nt]
    nt=np.where((np.abs(x_pos) > 10) & (np.abs(y_pos) > 10))
    fib_ids=fib_id[nt]
        
    #import matplotlib.pyplot as plt
    #plt.plot(x_posf*psc,y_posf*psc,'o')
    #plt.show()
    x_ifu=x_posf*psc
    y_ifu=y_posf*psc
    return x_ifu,y_ifu,fib_idt,fib_ids

def read_standar(path_data='',stdar_t='Feige32',stdT='',fergs=True):
    wav_i=[]
    res_i=[]
    file=path_data+'/'+stdar_t+stdT+'.dat'
    f=open(file,'r')
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=list(filter(None,data))
            wav_i.extend([float(data[0])])
            if fergs:
                res_i.extend([float(data[1])*1e-16])
            else:
                res_i.extend([float(data[1])])
    f.close()
    wav_i=np.array(wav_i)
    res_i=np.array(res_i)
    return wav_i,res_i