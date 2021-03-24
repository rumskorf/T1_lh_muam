# /usr/bin/python

'''                                                                    

HELP

Python Libraries for:

        DATA: reading UKMO,NCEP,SABER,COMMA,TEC ...
              generating date
        CALC: calculating Period-Wavenumber Spectra (PK)
              calculating EP-Flux
        PLOT: plotting running/averaged Spectra
        SAVE:

usage:

        import sys
        import sys.path.append('~/ANALYSIS/PYT/lib/')

        import cpw_lib as LIB
        
# -*- coding: utf-8 -*-
# f2py --fcompiler=gfortran -m -c grwaves grwaves.f
# 
# cd src
# grep c_p *.f
#

'''

__author__ = "Peter Hoffmann"
__date__ = "15 March 2009"
__version__ = "v1"
__credits__ = """Institute for Meteorology, University of Leipzig"""


#####################################################################

import numpy as N

import array
#import commands

from netCDF4 import Dataset

##############################################################################################

nx = 64
ny = 36
nz = 56
nv = 3 #anzahl variablen
nh = 24 #intervalle pro tag (6=4h, 12=2h)
nd = 30
nt = nh*nd
nw = nd
	
lon = N.arange(0.,5.625*nx,5.625,float)
lat = N.arange(-87.5,-87.5+5.*ny,5.,float)
lev = N.linspace(1.421,1.421+2.842*(nz-1),nz,float)
day = N.arange(0.,float(nw*nh),1,float)

print('lon',lon.shape)#,lon
print('lat',lat.shape)#,lat
print('lev',lev.shape)#,lev
print('day',day.shape)#,day

runs = 'Jan330'			# eintragen des file namens

fileobj=open('gwfluxu_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
gwu = N.array(binvalues, float)   
gwu = N.reshape(gwu,(nt,ny,nx,nz))
fileobj.close() 

print('gwu: ',gwu.shape)

fileobj=open('gwfluxv_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
gwv = N.array(binvalues, float)
gwv = N.reshape(gwv,(nt,ny,nx,nz))
fileobj.close() 

print('gwv: ',gwv.shape)

fileobj=open('gwacu_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
gau = N.array(binvalues, float)
gau = N.reshape(gau,(nt,ny,nx,nz))
fileobj.close() 

print('gau: ',gau.shape)

fileobj=open('gwacv_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
gav = N.array(binvalues, float)
gav = N.reshape(gav,(nt,ny,nx,nz))
fileobj.close() 

print('gav: ',gav.shape)

fileobj=open('gwt_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
gat = N.array(binvalues, float)
gat = N.reshape(gat,(nt,ny,nx,nz))
fileobj.close() 

print('gat: ',gat.shape)

  #########################################################################

t1 = 0*nh
t2 = nd*nh

gfu = N.zeros((nw*nh,nz,ny,nx),float)
gfv = N.zeros((nw*nh,nz,ny,nx),float)
gcu = N.zeros((nw*nh,nz,ny,nx),float)
gcv = N.zeros((nw*nh,nz,ny,nx),float)
ght = N.zeros((nw*nh,nz,ny,nx),float)

for z in range(nz):
    gfu[:,z,:,:] = gwu[t1:t2,:,:,z]
    gfv[:,z,:,:] = gwv[t1:t2,:,:,z]
    gcu[:,z,:,:] = gau[t1:t2,:,:,z]
    gcv[:,z,:,:] = gav[t1:t2,:,:,z]
    ght[:,z,:,:] = gat[t1:t2,:,:,z]

print('gfu: ',gfu.shape)
print('gfv: ',gfv.shape)
print('gcu: ',gcu.shape)
print('gcv: ',gcv.shape)
print('ght: ',ght.shape)

  #########################################################################

out = Dataset('muam_'+runs+'_gw.nc','w')
	  
out.createDimension('lons',nx)                
out.createDimension('lats',ny)
out.createDimension('levs',nz)
out.createDimension('time',nt)  
		    
v = out.createVariable('gfu','f',('time','levs','lats','lons'))
v[:,:,:,:]=gfu[:,:,:,:]
v.units = 'm2/s2'

v = out.createVariable('gfv','f',('time','levs','lats','lons'))
v[:,:,:,:]=gfv[:,:,:,:]
v.units = 'm2/s2'

v = out.createVariable('gcu','f',('time','levs','lats','lons'))
v[:,:,:,:]=gcu[:,:,:,:]
v.units = 'm/s/s'

v = out.createVariable('gcv','f',('time','levs','lats','lons'))
v[:,:,:,:]=gcv[:,:,:,:]
v.units = 'm/s/s'

v = out.createVariable('ght','f',('time','levs','lats','lons'))
v[:,:,:,:]=ght[:,:,:,:]
v.units = 'K/s'

v = out.createVariable('day','f',('time',))
v[:]=day[:]
v.units = 'h'

v = out.createVariable('lev','f',('levs',))
v[:]=lev[:]
v.units = 'km'

v = out.createVariable('lat','f',('lats',))
v[:]=lat[:]
v.units = 'deg'

v = out.createVariable('lon','f',('lons',))
v[:]=lon[:]
v.units = 'deg'

out.close()
###########
