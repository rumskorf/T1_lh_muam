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
nh = 24 #intervalle pro tag (6=4h, 12=2h, 24=1h)
nd = 30
nt = nh*nd
nw = nd
	
lon = N.arange(0.,5.625*nx,5.625,float)
lat = N.arange(-87.5,-87.5+5.*ny,5.,float)
lev = N.arange(1.421,160.,2.842,float)
day = N.arange(0.,float(nw*nh),1,float)

print('lon',lon.shape)#,lon
print('lat',lat.shape)#,lat
print('lev',lev.shape)#,lev
print('day',day.shape)#,day

runs = 'Jan330'			# eintragen des file namens

fileobj=open('uvt_'+runs+'.dx',mode='rb') 	
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt*nv)
uvt = N.array(binvalues, float)
uvt = N.reshape(uvt,(nt,nv,nz,ny,nx))
fileobj.close()

print('uvt: ',uvt.shape)

fileobj=open('phi_'+runs+'.dx',mode='rb')
binvalues = array.array('f') 
binvalues.fromfile(fileobj, nx*ny*nz*nt)
phi1 = N.array(binvalues, float)  
phi1 = N.reshape(phi1,(nt,nz,ny,nx))    
fileobj.close()

print('phi: ',phi1.shape)

fileobj=open('wvel_'+runs+'.dx',mode='rb')
binvalues = array.array('f')
binvalues.fromfile(fileobj, nx*ny*nz*nt)
ver1 = N.array(binvalues, float)
ver1 = N.reshape(ver1,(nt,nz,ny,nx))
fileobj.close()

print('ver: ',ver1.shape)

#########################################################################

t1 = 0*nh
t2 = nd*nh

zon = uvt[t1:t2,0,:,:,:]
mer = uvt[t1:t2,1,:,:,:]
tem = uvt[t1:t2,2,:,:,:]
phi = phi1[t1:t2,:,:,:]
ver = ver1[t1:t2,:,:,:]

#########################################################################

out = Dataset('muam_'+runs+'.nc','w')
                
out.createDimension('lons',len(lon))                
out.createDimension('lats',len(lat))
out.createDimension('levs',len(lev))
out.createDimension('time',len(day))
                          
v = out.createVariable('zon','f',('time','levs','lats','lons'))
v[:,:,:,:]=zon[:,:,:,:]
v.units = 'm/s'

v = out.createVariable('mer','f',('time','levs','lats','lons'))
v[:,:,:,:]=mer[:,:,:,:]
v.units = 'm/s'

v = out.createVariable('tem','f',('time','levs','lats','lons'))
v[:,:,:,:]=tem[:,:,:,:]
v.units = 'K'

v = out.createVariable('phi','f',('time','levs','lats','lons'))
v[:,:,:,:]=phi[:,:,:,:]
v.units = 'm'

v = out.createVariable('ver','f',('time','levs','lats','lons'))
v[:,:,:,:]=ver[:,:,:,:]
v.units = 'm/s'

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
