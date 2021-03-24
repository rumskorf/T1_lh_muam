#!/usr/bin/env python
# coding: utf-8

# In[33]:


import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pathlib import Path
import numpy as np


# In[72]:


xr.set_options(keep_attrs = True) 
plt.rcParams.update({'font.size': 20})

c_pad = 0.2
c_shrink = 0.99
c_aspect = 20

ymajorLocator = MultipleLocator(20)
xmajorLocator = MultipleLocator(30)


# In[96]:


def process(ds):
    ds['lat'] = lat_grid
    ds['lat'].attrs['long_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'deg'
    ds['lev'] = lev_grid
    ds['lev'].attrs['long_name'] = 'height'
    ds['lev'].attrs['units'] = 'km'
    return ds.mean('yr'), ds.std('yr')

def phase_mean(x,y):
    xsum = x.sum('yr')
    ysum = y.sum('yr')
    ampsum=np.sqrt(xsum**2+ysum**2)
    phamean=np.arctan2(ysum,xsum)
    return phamean

def process_vec(ds, sel_tide):
    ds['lat'] = lat_grid
    ds['lat'].attrs['long_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'deg'
    ds['lev'] = lev_grid
    ds['lev'].attrs['long_name'] = 'height'
    ds['lev'].attrs['units'] = 'km'
    
    x = ds[sel_tide+'_amp']*np.cos(ds[sel_tide+'_pha'])
    y = ds[sel_tide+'_amp']*np.sin(ds[sel_tide+'_pha'])

    return phase_mean(x,y), x, y


# In[3]:


root_path = '/home/gemeinsam_tmp/VACILT/'

ny = 36
nz = 56
nt = 5
lev_grid = np.arange(1.421,160.,2.842,float)
lat_grid = np.arange(-87.5,-87.5+5.*ny,5.,float)


# In[59]:

sel_var_array = ['tem', 'zon', 'mer', 'ver', 'phi']
# sel_var = 'tem'
for sel_var in sel_var_array:
    pathh = str('/home/gemeinsam_tmp/VACILT/latent_heat_output/comparison2/'+str(sel_var))
    sel_month = 'Jan'
    sel_phase_array = ['la','el']
    # sel_phase = 'el'
    #for sel_phase in sel_phase_array:


    # In[60]:


#    infile = f'ensemble_tides_{sel_var}_{sel_phase}_{sel_month}.nc'
    
    infile_wo_el = f'ensemble_tides_{sel_var}_el_{sel_month}.nc'
    infile_wo_la = f'ensemble_tides_{sel_var}_la_{sel_month}.nc'
    infile_w_el = f'ensemble_tides_{sel_var}_el_{sel_month}.nc'
    infile_w_la = f'ensemble_tides_{sel_var}_la_{sel_month}.nc'
    
    #ds_wo_lh = xr.open_dataset(f'{root_path}muam_mstober/{infile}')
    #ds_w_lh = xr.open_dataset(f'{root_path}latent_heat_output/{infile}')
    
    el_wo_lh = xr.open_dataset(f'{root_path}muam_mstober/{infile_wo_el}')
    la_wo_lh = xr.open_dataset(f'{root_path}muam_mstober/{infile_wo_la}')
    el_w_lh = xr.open_dataset(f'{root_path}muam_mstober/{infile_w_el}')
    la_w_lh = xr.open_dataset(f'{root_path}muam_mstober/{infile_w_la}')


    # In[121]:


    # sel_tide = 'DT'
    sel_tide_array = ['DT','SDT','TDT','QDT']
    for sel_tide in sel_tide_array:
        # sel_ave = 'amp'
        sel_ave_array = ['amp','pha']
        for sel_ave in sel_ave_array:


            # In[122]:


            if sel_ave == 'amp':
                #mean_wo_lh, std_wo_lh = process(ds_wo_lh[f'{sel_tide}_{sel_ave}'])
                #mean_w_lh, std_w_lh = process(ds_w_lh[f'{sel_tide}_{sel_ave}'])
                mean_wo_lh_el, std_wo_lh_el = process(el_wo_lh[f'{sel_tide}_{sel_ave}'])
                mean_wo_lh_la, std_wo_lh_la = process(la_wo_lh[f'{sel_tide}_{sel_ave}'])
                mean_w_lh_el, std_w_lh_el = process(el_w_lh[f'{sel_tide}_{sel_ave}'])
                mean_w_lh_la, std_w_lh_la = process(la_w_lh[f'{sel_tide}_{sel_ave}'])
                #diff = mean_w_lh-mean_wo_lh
                diff_wo = mean_wo_lh_el-mean_wo_lh_la
                diff_w = mean_w_lh_el-mean_w_lh_la
                diff_diff = diff_w-diff_wo
                
                std_wo_lh = std_wo_lh_el + std_wo_lh_la
                std_w_lh = std_w_lh_el + std_w_lh_la

                cmap=plt.get_cmap('Blues') #cmap with blue colors
                #make descrete areas in colorbar
                cmaplist=[cmap(i) for i in range(cmap.N)]
                cmaplist[0]='white'
                cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N) 
                
                if sel_var=='tem':
                    print("tem")
                    dsig=0.4
                    bounds = [0.0, 0.5, 1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0] # colorbar range bounds = N.arange(0,30.1,2)
                    lvls = np.arange(dsig,5,dsig)
                    diff_clabel = 'diff [K]'
                    clabel = 'temperature [K]'
                elif sel_var=='zon':
                    print("zon")
                    dsig=0.5
                    bounds = np.arange(0,36.1,3) # colorbar range
                    lvls = np.arange(dsig,5,dsig)
                    diff_clabel = 'diff [m/s]'
                    clabel = 'zonal wind [m/s]'
                elif sel_var=='mer':
                    print("mer")
                    dsig=0.5
                    bounds = np.arange(0,36.1,3) # colorbar range
                    lvls = np.arange(dsig,5,dsig)
                    diff_clabel = 'diff [m/s]'
                    clabel = 'meridional wind [m/s]'
                elif sel_var=='ver':
                    print("ver")
                    dsig=0.2
                    bounds = np.arange(0,10.1,1) # colorbar range
                    lvls = np.arange(dsig,5,dsig)
                    diff_clabel = 'diff [m/s]'
                    clabel = 'vertical wind [m/s]'
                elif sel_var=='phi':
                    print("phi")
                    dsig=5
                    bounds = np.arange(0,1001,100) # colorbar range
                    lvls = np.arange(dsig,50,dsig)
                    diff_clabel = 'diff [m]'
                    clabel = 'geopotential height [m]'
                    
                cmap.set_over('black')
                

                
            elif sel_ave == 'pha':
                #mean_wo_lh, x_wo_lh, y_wo_lh = process_vec(ds_wo_lh, sel_tide)
                #mean_w_lh, x_w_lh, y_w_lh = process_vec(ds_w_lh, sel_tide)
                mean_wo_lh_el, x_wo_lh_el, y_wo_lh_el = process_vec(el_wo_lh, sel_tide)
                mean_wo_lh_la, x_wo_lh_la, y_wo_lh_la = process_vec(la_wo_lh, sel_tide)
                mean_w_lh_el, x_w_lh_el, y_w_lh_el = process_vec(el_w_lh, sel_tide)
                mean_w_lh_la, x_w_lh_la, y_w_lh_la = process_vec(la_w_lh, sel_tide)
                #diff_x = x_w_lh-x_wo_lh
                #diff_y = y_w_lh-y_wo_lh
                #diff = mean_w_lh-mean_wo_lh #phase_mean(diff_x, diff_y)   
                diff_wo = mean_wo_lh_el-mean_wo_lh_la #phase_mean(diff_x, diff_y)   
                diff_w = mean_w_lh_el-mean_w_lh_la #phase_mean(diff_x, diff_y)   
                diff_diff = diff_w-diff_wo #phase_mean(diff_x, diff_y)   
                
                cmap=plt.get_cmap('Reds')
                cmaplist=[cmap(i) for i in range(cmap.N)]
                cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N) 
                bounds = np.arange(-np.pi,1.1*np.pi,np.pi/2.) # colorbar range
                clabel = 'phase [rad]'
                diff_clabel = 'diff [rad]'  


            # In[130]:


            plt.rcParams.update({'font.size': 20})
            fig, axes = plt.subplots(ncols = 3, figsize = (16,6), sharey = True, sharex = True)

            ax = axes[0]
            p = diff_w.plot(cmap = cmap, ax = ax, 
                                         add_colorbar = False) #, levels = bounds
            if sel_ave == 'amp':
                std_w_lh.squeeze().plot.contour(colors='g', ax = ax) #levels = lvls,
                ax.text(0.6,0.1,'$\sigma_{max}$'+f'= {std_w_lh.max().values:2.3f}', 
                    fontsize=13, backgroundcolor=(1, 1, 1, 0.5), transform=ax.transAxes)


            ax.set_title('Diff El-La with LH')
            ax.xaxis.set_major_locator(xmajorLocator)
            ax.yaxis.set_major_locator(ymajorLocator)
            ax.grid(b=True, which='major', color='gray', linestyle='--')

            ax = axes[1]
            cbar = fig.colorbar(p, ax=axes[:2],  location='bottom', extend = 'both', pad = c_pad,  shrink=c_shrink, aspect = c_aspect)
            cbar.set_label(clabel)
            #if sel_ave == 'pha':
                #cbar.ax.set_xticklabels((r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$',r'$\pi$')) 

            diff_wo.plot(cmap = cmap, ax = ax,
                                      add_colorbar = False) #, levels = bounds
            if sel_ave == 'amp':
                std_wo_lh.squeeze().plot.contour(colors='g', ax = ax) #levels = lvls,
                ax.text(0.6,0.1,'$\sigma_{max}$'+f'= {std_wo_lh.max().values:2.3f}', 
                    fontsize=13, backgroundcolor=(1, 1, 1, 0.5), transform=ax.transAxes)

            ax.set_ylabel('')
            ax.set_title('Diff El-La without LH')
            ax.grid(b=True, which='major', color='gray', linestyle='--')

            ax = axes[2]
            p = diff_diff.plot(ax = ax, robust = True, 
                                cbar_kwargs={'orientation': 'horizontal', 'pad': c_pad, \
                                             'aspect': c_aspect/2., "shrink": c_shrink, 
                                             'extend': 'both', 'label': diff_clabel})
            ax.text(0.6,0.1,'$\Delta_{max}$'+f'= {diff_diff.max().values:2.3f}', 
                    fontsize=13, backgroundcolor=(1, 1, 1, 0.5), transform=ax.transAxes)
            #if sel_ave == 'pha':
                #p.colorbar.set_ticks(bounds/10.)
                #p.colorbar.set_ticklabels((r'$-\pi/10$',r'$-\pi/20$', '0', r'$\pi/20$',r'$\pi/10$'))
                
            ax.set_ylabel('')
            ax.set_title('Diff of Diff (with-without)')
            ax.grid(b=True, which='major', color='gray', linestyle='--')

            plt.suptitle(f'{sel_tide} comparison for {sel_month}', y=1.02)
            fig.savefig(f'{pathh}/{sel_tide}_{sel_ave}_comparison2_{sel_var}_{sel_month}.png', bbox_inches='tight')
                        
