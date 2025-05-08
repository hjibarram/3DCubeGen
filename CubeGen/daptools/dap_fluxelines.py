from astropy.io import fits
import numpy as np
from tqdm.notebook import tqdm
from tqdm import tqdm as tqdmT
from astropy.wcs import WCS
import CubeGen.tools.tools as tools

def extract_values(list_vals,hdu_list,key_list='NP_ELINES_B'):
    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention=list_vals
    table_hdu=hdu_list[key_list]
    table_hdr=table_hdu.header
    keys=list(table_hdr.keys())
    for it in range(0, len(keys)):
        if 'TTYPE' in keys[it]:
            extention.extend([key_list])
            if 'flux' in table_hdr[keys[it]] and not 'e_flux' in table_hdr[keys[it]]:
                keys_flux.extend([table_hdr[keys[it]]])
            if 'e_flux' in table_hdr[keys[it]]:
                keys_eflux.extend([table_hdr[keys[it]]])
            if 'disp' in table_hdr[keys[it]] and not 'e_disp' in table_hdr[keys[it]]:
                keys_disp.extend([table_hdr[keys[it]]])
            if 'e_disp' in table_hdr[keys[it]]:
                keys_edisp.extend([table_hdr[keys[it]]])
            if 'vel' in table_hdr[keys[it]] and not 'e_vel' in table_hdr[keys[it]]:
                keys_vel.extend([table_hdr[keys[it]]])
            if 'e_vel' in table_hdr[keys[it]]:
                keys_evel.extend([table_hdr[keys[it]]])
            if 'EW' in table_hdr[keys[it]] and not 'e_EW' in table_hdr[keys[it]]:
                keys_EW.extend([table_hdr[keys[it]]])
            if 'e_EW' in table_hdr[keys[it]]:
                keys_eEW.extend([table_hdr[keys[it]]])    
    return [keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention]


def extract_PMvalues(list_vals,hdu_list,key_list='PM_ELINES'):
    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,extention=list_vals
    table_hdu=hdu_list[key_list]
    table_hdr=table_hdu.header
    keys=list(table_hdr.keys())
    for it in range(0, len(keys)):
        if 'TTYPE' in keys[it]:
            extention.extend([key_list])
            if 'flux' in table_hdr[keys[it]] and not 'e_flux' in table_hdr[keys[it]]:
                keys_flux.extend([table_hdr[keys[it]]])
            if 'e_flux' in table_hdr[keys[it]]:
                keys_eflux.extend([table_hdr[keys[it]]])
            if 'disp' in table_hdr[keys[it]] and not 'e_disp' in table_hdr[keys[it]]:
                keys_disp.extend([table_hdr[keys[it]]])
            if 'e_disp' in table_hdr[keys[it]]:
                keys_edisp.extend([table_hdr[keys[it]]])
            if 'vel' in table_hdr[keys[it]] and not 'e_vel' in table_hdr[keys[it]]:
                keys_vel.extend([table_hdr[keys[it]]])
            if 'e_vel' in table_hdr[keys[it]]:
                keys_evel.extend([table_hdr[keys[it]]])
    return [keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,extention]

def get_values(file_rss,file_dap,nx,ny,map_dap=[],map_dap2=[],notebook=True):
    
    hdu_list = fits.open(file_rss)
    table_hdu = hdu_list['SLITMAP']
    table_data = table_hdu.data
    ypt0=table_data.field('xpmm')
    xpt0=table_data.field('ypmm')
    ra_fib=table_data.field('ra')
    dec_fib=table_data.field('dec')
    Std_id=table_data.field('fiberid')-1
    hdu_list.close()
    
    hdu_list = fits.open(file_dap)
    extention=[]
    keys_disp=[]
    keys_vel=[]
    keys_flux=[]
    keys_EW=[]
    keys_edisp=[]
    keys_evel=[]
    keys_eflux=[]
    keys_eEW=[]
    
    list_vals=[keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention]
    list_vals_out=extract_values(list_vals,hdu_list,key_list='NP_ELINES_B')
    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention=list_vals_out
    list_vals=[keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention]
    list_vals_out=extract_values(list_vals,hdu_list,key_list='NP_ELINES_R')
    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention=list_vals_out
    list_vals=[keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention]
    list_vals_out=extract_values(list_vals,hdu_list,key_list='NP_ELINES_I')
    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW,extention=list_vals_out
    n_keys=len(keys_flux)
    
    table_hduB = hdu_list['NP_ELINES_B']
    table_dataB = table_hduB.data
    table_hduR = hdu_list['NP_ELINES_R']
    table_dataR = table_hduR.data
    table_hduI = hdu_list['NP_ELINES_I']
    table_dataI = table_hduI.data

    extention1=[]
    keys_disp1=[]
    keys_vel1=[]
    keys_flux1=[]
    keys_edisp1=[]
    keys_evel1=[]
    keys_eflux1=[]
    list_vals1=[keys_flux1,keys_eflux1,keys_disp1,keys_edisp1,keys_vel1,keys_evel1,extention1]
    list_vals1_out=extract_PMvalues(list_vals1,hdu_list,key_list='PM_ELINES')
    keys_flux1,keys_eflux1,keys_disp1,keys_edisp1,keys_vel1,keys_evel1,extention1=list_vals1_out
    n_keys1=6
    keys_vals=[keys_flux1[0],keys_vel1[0],keys_disp1[0],keys_eflux1[0],keys_evel1[0],keys_edisp1[0]]
    table_hduPM = hdu_list['PM_ELINES']
    table_dataPM = table_hduPM.data

    id_PM=table_dataPM.field('id_fib')
    wave_mod=table_dataPM.field('wl')
    nt=np.where(id_PM == 0)[0]
    n_mods=len(nt)
    wave_mod=wave_mod[nt]
    id_dap=table_dataB.field('id')    
    xp=[]
    yp=[]
    if notebook:
        pbar=tqdm(total=len(id_dap)*len(Std_id))
    else:
    	pbar=tqdmT(total=len(id_dap)*len(Std_id))
    for j in range(0, len(id_dap)):
        for k in range(0, len(Std_id)):
            if int(id_dap[j].split('.')[1])-1 == Std_id[k]:
                xp.extend([xpt0[k]])
                yp.extend([ypt0[k]])
            pbar.update(1)
    pbar.close()
    xp=np.array(xp)
    yp=np.array(yp)

    n_spax=len(xp)
    flux_val=np.zeros([n_keys,n_spax])
    disp_val=np.zeros([n_keys,n_spax])
    vel_val=np.zeros([n_keys,n_spax])
    EW_val=np.zeros([n_keys,n_spax])
    eflux_val=np.zeros([n_keys,n_spax])
    edisp_val=np.zeros([n_keys,n_spax])
    evel_val=np.zeros([n_keys,n_spax])
    eEW_val=np.zeros([n_keys,n_spax])

    flux_valPM=np.zeros([n_mods,n_spax])
    disp_valPM=np.zeros([n_mods,n_spax])
    vel_valPM=np.zeros([n_mods,n_spax])
    eflux_valPM=np.zeros([n_mods,n_spax])
    edisp_valPM=np.zeros([n_mods,n_spax])
    evel_valPM=np.zeros([n_mods,n_spax])

    
    for j in range(0, n_keys):
        try:
            val1=table_dataB.field(keys_flux[j])
            val2=table_dataB.field(keys_disp[j])
            val3=table_dataB.field(keys_vel[j])
            val4=table_dataB.field(keys_EW[j])
            val1e=table_dataB.field(keys_eflux[j])
            val2e=table_dataB.field(keys_edisp[j])
            val3e=table_dataB.field(keys_evel[j])
            val4e=table_dataB.field(keys_eEW[j])
        except:
            try:
                val1=table_dataR.field(keys_flux[j])
                val2=table_dataR.field(keys_disp[j])
                val3=table_dataR.field(keys_vel[j])
                val4=table_dataR.field(keys_EW[j])
                val1e=table_dataR.field(keys_eflux[j])
                val2e=table_dataR.field(keys_edisp[j])
                val3e=table_dataR.field(keys_evel[j])
                val4e=table_dataR.field(keys_eEW[j])
            except:
                val1=table_dataI.field(keys_flux[j])
                val2=table_dataI.field(keys_disp[j])
                val3=table_dataI.field(keys_vel[j])
                val4=table_dataI.field(keys_EW[j])
                val1e=table_dataI.field(keys_eflux[j])
                val2e=table_dataI.field(keys_edisp[j])
                val3e=table_dataI.field(keys_evel[j])
                val4e=table_dataI.field(keys_eEW[j])
        flux_val[j,:]=val1
        disp_val[j,:]=val2
        vel_val[j,:]=val3
        EW_val[j,:]=val4
        eflux_val[j,:]=val1e
        edisp_val[j,:]=val2e
        evel_val[j,:]=val3e
        eEW_val[j,:]=val4e
    for j in range(0, n_spax):
        nt=np.where(id_PM == j)[0]
        val1=table_dataPM.field(keys_flux1[0])[nt]
        val2=table_dataPM.field(keys_disp1[0])[nt]
        val3=table_dataPM.field(keys_vel1[0])[nt]
        val1e=table_dataPM.field(keys_eflux1[0])[nt]
        val2e=table_dataPM.field(keys_edisp1[0])[nt]
        val3e=table_dataPM.field(keys_evel1[0])[nt]
        flux_valPM[:,j]=val1
        disp_valPM[:,j]=val2
        vel_valPM[:,j]=val3
        eflux_valPM[:,j]=val1e
        edisp_valPM[:,j]=val2e
        evel_valPM[:,j]=val3e
    hdu_list.close()

    if len(map_dap) == 0:
        map_dap=np.zeros([n_keys*8,nx,ny])
    if len(map_dap2) == 0:
        map_dap2=np.zeros([n_mods*6,nx,ny])    
    
    if notebook:    
        pbar=tqdm(total=n_spax)
    else:
    	pbar=tqdmT(total=n_spax)
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, n_spax):
                if i==xp[k] and j == yp[k]:
                    for h in range(0, n_keys):
                        map_dap[h+n_keys*0,i,j]=flux_val[h,k]
                        map_dap[h+n_keys*1,i,j]=vel_val[h,k]
                        map_dap[h+n_keys*2,i,j]=disp_val[h,k]
                        map_dap[h+n_keys*3,i,j]=EW_val[h,k]
                        map_dap[h+n_keys*4,i,j]=eflux_val[h,k]
                        map_dap[h+n_keys*5,i,j]=evel_val[h,k]
                        map_dap[h+n_keys*6,i,j]=edisp_val[h,k]
                        map_dap[h+n_keys*7,i,j]=eEW_val[h,k]
                    for h in range(0, n_mods):
                        map_dap2[h+n_mods*0,i,j]=flux_valPM[h,k]
                        map_dap2[h+n_mods*1,i,j]=vel_valPM[h,k]
                        map_dap2[h+n_mods*2,i,j]=disp_valPM[h,k]
                        map_dap2[h+n_mods*3,i,j]=eflux_valPM[h,k]
                        map_dap2[h+n_mods*4,i,j]=evel_valPM[h,k]
                        map_dap2[h+n_mods*5,i,j]=edisp_valPM[h,k]
                    pbar.update(1)
    pbar.close()

    list_vals=[keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW]
    
    return map_dap,map_dap2,n_keys,list_vals,wave_mod,keys_vals

def dap_extract_fluxelines(name,out_path='./',path_cube='./',path_dap='./',path_rss='./',basename_cube='lvmCube-NAME.fits.gz',basename_dap='NAME.dap.fits.gz',basename_rss='lvmRSS-NAME.fits.gz',nsplit=0,spt=[False,[1,1],2],notebook=True):    

    if nsplit > 1:
        file=path_cube+'/'+basename_cube.replace('NAME',name+'_00')
        hdr=fits.getheader(file, 0)
    else:
        file=path_cube+'/'+basename_cube.replace('NAME',name)
        hdr=fits.getheader(file, 0)
    ny=hdr['NAXIS1']
    nx=hdr['NAXIS2']
    wcs=WCS(hdr).celestial
    head = wcs.to_header()
    
    
    if nsplit > 1:
        for i in range(0, nsplit):
            for j in range(0, nsplit):
                val='_p'+str(i)+str(j)
                file_rss=path_rss+'/'+basename_rss.replace('NAME',name+val)
                file_dap=path_dap.replace(name,name+val)+'/'+basename_dap.replace('NAME',name+val)
                if i == 0 and j == 0:
                    map_dap,map_dap2,n_keys,list_vals_out,wave_mod,keys_vals=get_values(file_rss,file_dap,nx,ny,notebook=notebook)
                else:
                    #if i == 2 and j ==3:
                    if spt[0] == True:
                        if spt[1][0] == i and spt[1][1] == j:
                            for k in range(0, spt[2]):
                                for l in range(0, spt[2]):
                                    val='_p'+str(i)+str(j)+'_p'+str(k)+str(l)
                                    file_rss=path_rss+'/'+basename_rss.replace('NAME',name+val)
                                    file_dap=path_dap.replace(name,name+val)+'/'+basename_dap.replace('NAME',name+val)
                                    map_dap,map_dap2,n_keys,list_vals_out=get_values(file_rss,file_dap,nx,ny,map_dap=map_dap,map_dap2=map_dap2,notebook=notebook)
                    else:
                        map_dap,map_dap2,n_keys,list_vals_out,wave_mod,keys_vals=get_values(file_rss,file_dap,nx,ny,map_dap=map_dap,map_dap2=map_dap2,notebook=notebook)
    else:
        file_rss=path_rss+'/'+basename_rss.replace('NAME',name)
        file_dap=path_dap+'/'+basename_dap.replace('NAME',name)
        map_dap,map_dap2,n_keys,list_vals_out,wave_mod,keys_vals=get_values(file_rss,file_dap,nx,ny,notebook=notebook)

    keys_flux,keys_eflux,keys_disp,keys_edisp,keys_vel,keys_evel,keys_EW,keys_eEW=list_vals_out
    
    tools.sycall('mkdir -p '+out_path) 

    h1=fits.PrimaryHDU(map_dap)#,header=head)
    h2=fits.ImageHDU(map_dap2)
    h_k=h1.header
    h_k["CRVAL1"]=hdr['CRVAL1']
    h_k["CDELT1"]=hdr["CD1_1"]
    h_k["CD1_1"]=hdr["CD1_1"]
    h_k["CD1_2"]=hdr["CD1_2"]
    h_k["CRPIX1"]=hdr["CRPIX1"]
    h_k["CTYPE1"]=hdr["CTYPE1"]
    h_k["CRVAL2"]=hdr['CRVAL2']
    h_k["CD2_1"]=hdr["CD2_1"]
    h_k["CD2_2"]=hdr["CD2_2"]
    h_k["CDELT2"]=hdr["CD2_2"]
    h_k["CRPIX2"]=hdr["CRPIX2"]
    h_k["CTYPE2"]=hdr["CTYPE2"]
    h_k["CUNIT1"]=hdr["CUNIT1"]                                           
    h_k["CUNIT2"]=hdr["CUNIT2"]
    h_k["RADESYS"]=hdr["RADESYS"]
    h_k["OBJSYS"]=hdr["OBJSYS"]
    h_k["EQUINOX"]=hdr["EQUINOX"]    
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*0)))
        h_k["VAL_"+lab]=keys_flux[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*1)))
        h_k["VAL_"+lab]=keys_vel[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*2)))
        h_k["VAL_"+lab]=keys_disp[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*3)))
        h_k["VAL_"+lab]=keys_EW[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*4)))
        h_k["VAL_"+lab]=keys_eflux[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*5)))
        h_k["VAL_"+lab]=keys_evel[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*6)))
        h_k["VAL_"+lab]=keys_edisp[i]
    for i in range(0, n_keys):
        lab="{:0>4}".format(str(int(i+1+n_keys*7)))
        h_k["VAL_"+lab]=keys_eEW[i]
    h_k['EXTNAME']='NP'    
    h_k.update()
    n_mods=len(wave_mod)
    n_keys2=len(keys_vals)
    h_k=h2.header
    h_k["CRVAL1"]=hdr['CRVAL1']
    h_k["CD1_1"]=hdr["CD1_1"]
    h_k["CD1_2"]=hdr["CD1_2"]
    h_k["CRPIX1"]=hdr["CRPIX1"]
    h_k["CTYPE1"]=hdr["CTYPE1"]
    h_k["CRVAL2"]=hdr['CRVAL2']
    h_k["CD2_1"]=hdr["CD2_1"]
    h_k["CD2_2"]=hdr["CD2_2"]
    h_k["CRPIX2"]=hdr["CRPIX2"]
    h_k["CTYPE2"]=hdr["CTYPE2"]
    h_k["CUNIT1"]=hdr["CUNIT1"]                                           
    h_k["CUNIT2"]=hdr["CUNIT2"]
    h_k["RADESYS"]=hdr["RADESYS"]
    h_k["OBJSYS"]=hdr["OBJSYS"]
    h_k["EQUINOX"]=hdr["EQUINOX"]  
    h_k['EXTNAME']='PM'
    for j in range(0, n_keys2):
        for i in range(0, n_mods):
            lab="{:0>4}".format(str(int(i+1+n_mods*j)))
            h_k["VAL_"+lab]=str(np.round(wave_mod[i],2))+'_'+keys_vals[j]
    h_k.update()
    hlist=fits.HDUList([h1,h2])
    hlist.update_extend()
    basenameC='lvmDAPMap-NAME-flux_elines.fits'
    file=out_path+'/'+basenameC.replace('NAME',name)
    out_fit=file
    hlist.writeto(out_fit,overwrite=True)
    tools.sycall('gzip -f '+out_fit)