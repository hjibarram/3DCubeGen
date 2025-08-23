import numpy as np

def read_vpwcs(name,path_data='data',base_name='NAME_wcs.txt',hdrid='#'):
    fibid=[]
    ra=[]
    dec=[]
    file=path_data+base_name.replace('NAME',name)
    f=open(file,'r')
    for line in f:
        if not hdrid in line:
            data=line.replace('\n','').split(' ')
            data=list(filter(None,data))
            fibid.extend([int(data[0])])
            ra.extend([float(data[1])*3600.0])
            dec.extend([float(data[2])*3600.0])
    f.close()
    fibid=np.array(fibid)
    ra=np.array(ra)
    dec=np.array(dec)
    return fibid,ra,dec