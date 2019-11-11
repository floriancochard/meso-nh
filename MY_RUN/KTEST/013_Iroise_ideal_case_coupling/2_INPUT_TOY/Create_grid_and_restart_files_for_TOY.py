#!/usr/bin/python
# -*- coding: utf-8 -*-
#
###################################################
#=================================================#
#  Creating grid and restart file for toy model
#  Author : J. Pianezze
#  Date   :        2015
#=================================================#
###################################################
#
import netCDF4
import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
from pylab import *
import os
#
curdir_path=os.path.abspath(os.curdir)+'/'
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# To be defined by the user
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#-- Limit of the grid (from etopo)
lat_domain=[47.0, 49.5]
lon_domain=[-6.2, -4.0]

#-- Type of forcing to create the restart file 
#   for the toy : CNSTE or SINUS

# CNSTE
value_CNSTE=290.0

# SINUS
value_SINUS_COEF=0.011
value_SINUS_LENGTH=1000.
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

file_topo = netCDF4.Dataset('topo.nc')

#------ Read variables
lon_full=file_topo.variables['lon'][:]
lat_full=file_topo.variables['lat'][:]
topo_full=file_topo.variables['topo'][:,:]

ind_min_lon=find(lon_full[:]>lon_domain[0])[0] ; print ind_min_lon
ind_max_lon=find(lon_full[:]>lon_domain[1])[0] ; print ind_max_lon
ind_min_lat=find(lat_full[:]>lat_domain[0])[0] ; print ind_min_lat
ind_max_lat=find(lat_full[:]>lat_domain[1])[0] ; print ind_max_lat

lon=lon_full[ind_min_lon:ind_max_lon]
lat=lat_full[ind_min_lat:ind_max_lat]
topo=topo_full[ind_min_lat:ind_max_lat,ind_min_lon:ind_max_lon]

#plt.contourf(topo)
#plt.show()

nlon=np.size(lon) ;  print 'nlon=', nlon
nlat=np.size(lat) ;  print 'nlat=', nlat
ncorn=4           ;  print 'ncorn=', ncorn

print '---- longitude/latitude'
lon2D=np.zeros((nlat,nlon))
lat2D=np.zeros((nlat,nlon))

for ind_lon in xrange(nlon):
  lat2D[:,ind_lon]=lat[:]
for ind_lat in xrange(nlat):
  lon2D[ind_lat,:]=lon[:]

print '---- corners longitude/latitude'
clo=np.zeros((ncorn,nlat,nlon))
cla=np.zeros((ncorn,nlat,nlon))

deltax=lon[1]-lon[0] ; print 'deltax=', deltax
clo[0,:,:]=lon2D[:,:]+deltax/2.0
clo[1,:,:]=lon2D[:,:]-deltax/2.0
clo[2,:,:]=lon2D[:,:]-deltax/2.0
clo[3,:,:]=lon2D[:,:]+deltax/2.0

deltay=lat[1]-lat[0] ; print 'deltay=', deltay
cla[0,:,:]=lat2D[:,:]+deltay/2.0
cla[1,:,:]=lat2D[:,:]+deltay/2.0
cla[2,:,:]=lat2D[:,:]-deltay/2.0
cla[3,:,:]=lat2D[:,:]-deltay/2.0

print '---- surface'
surface=np.zeros((nlat,nlon))
surface[:,:]=deltax*deltay


print '---- mask and var send by toy'
mask=np.zeros((nlat,nlon))
toyvarcnste=np.zeros((nlat,nlon))
toyvarsinus=np.zeros((nlat,nlon))

for ind_lon in xrange(nlon):
  for ind_lat in xrange(nlat):
    if topo[ind_lat,ind_lon] > 0.0 :
      mask[ind_lat,ind_lon]=0
      toyvarcnste[ind_lat,ind_lon] = value_CNSTE
      toyvarsinus[ind_lat,ind_lon] = value_SINUS_COEF*math.sin(lat[ind_lat]*math.pi/180.0*value_SINUS_LENGTH)
    else:
      mask[ind_lat,ind_lon]=1
      toyvarcnste[ind_lat,ind_lon]=value_CNSTE
      toyvarsinus[ind_lat,ind_lon]= value_SINUS_COEF*math.sin(lat[ind_lat]*math.pi/180.0*value_SINUS_LENGTH)


##################################################
print '------------------------------------------'
print ' Creating netcdf file : grid_toy_model.nc'

grid_file=netCDF4.Dataset(curdir_path+'grid_toy_model.nc','w',format='NETCDF3_64BIT')
grid_file.Description='Grid file for OASIS coupling'

# ----------------------------------
# Create the dimensions of the files
# ----------------------------------
grid_file.createDimension ('nlon', nlon)
grid_file.createDimension ('nlat', nlat)
grid_file.createDimension ('ncorner', 4 )

# ----------------------------------
# Create the variables of the files
# ----------------------------------
varout=grid_file.createVariable('lon','d',('nlat','nlon'))
varout=grid_file.createVariable('lat','d',('nlat','nlon'))
varout=grid_file.createVariable('clo','d',('ncorner','nlat','nlon'))
varout=grid_file.createVariable('cla','d',('ncorner','nlat','nlon'))
varout=grid_file.createVariable('srf','d',('nlat','nlon'))
varout=grid_file.createVariable('imask','d',('nlat','nlon'))

# ---------------------------------------
# Write out the data arrays into the file
# ---------------------------------------
grid_file.variables['lon'][:,:] = lon2D[:,:]
grid_file.variables['lat'][:,:] = lat2D[:,:]
grid_file.variables['clo'][:,:] = clo[:,:,:]
grid_file.variables['cla'][:,:] = cla[:,:,:]
grid_file.variables['srf'][:,:] = surface[:,:]
grid_file.variables['imask'][:,:] = mask[:,:]

# ---------------------------------------
# close the file
# ---------------------------------------
grid_file.close()

print ' Closing netcdf file : grid_toy_model.nc'
print '------------------------------------------'
##################################################

##################################################
print '------------------------------------------'
print ' Creating netcdf file : rstrt_TOY.nc'

rstrt_file=netCDF4.Dataset(curdir_path+'rstrt_TOY.nc','w',format='NETCDF3_64BIT')
rstrt_file.Description='Restart file for OASIS coupling'

# ----------------------------------
# Create the dimensions of the files
# ----------------------------------
rstrt_file.createDimension ('nlon', nlon)
rstrt_file.createDimension ('nlat', nlat)

# ----------------------------------
# Create the variables of the files
# ----------------------------------
varout=rstrt_file.createVariable('VARCNSTE','d',('nlat','nlon'))
varout=rstrt_file.createVariable('VARSIN01','d',('nlat','nlon'))
varout=rstrt_file.createVariable('VARSIN02','d',('nlat','nlon'))

# ---------------------------------------
# Write out the data arrays into the file
# ---------------------------------------
rstrt_file.variables['VARCNSTE'][:,:] = toyvarcnste[:,:]
rstrt_file.variables['VARSIN01'][:,:] = toyvarsinus[:,:]
rstrt_file.variables['VARSIN02'][:,:] = toyvarsinus[:,:]

# ---------------------------------------
# close the file
# ---------------------------------------
rstrt_file.close()

print ' Closing netcdf file : rstrt_TOY.nc'
print '-----------------------------------------'
#####################################################
