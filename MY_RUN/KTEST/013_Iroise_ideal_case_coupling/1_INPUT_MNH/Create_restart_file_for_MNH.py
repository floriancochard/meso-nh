#!/usr/bin/python
# -*- coding: utf-8 -*-
#
##############################################
#============================================#
#      Creating restart file from prep real
#                 case files : mesonh model
#      Author : J. Pianezze
#      Date   :        2015
#============================================#
##############################################

import netCDF4
import numpy as np
import os

curdir_path=os.path.abspath(os.curdir)+'/'

file_RSTRT = netCDF4.Dataset(curdir_path+'IROISE_5KM_201109_02_00.nc')

LON_MNH=file_RSTRT.variables['LON'][1:-1,1:-1]
LAT_MNH=file_RSTRT.variables['LAT'][1:-1,1:-1]
U10_MNH=file_RSTRT.variables['UT'][2,1:-1,1:-1]
V10_MNH=file_RSTRT.variables['VT'][2,1:-1,1:-1]
PRES_MNH=file_RSTRT.variables['PABST'][2,1:-1,1:-1]

try:
  EVAP_MNH=file_RSTRT.variables['EVAP3D'][1:-1,1:-1]
except KeyError:
  print 'EVAP3D not found... imposed at 0!'
  EVAP_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))

try:
  RAIN_MNH=file_RSTRT.variables['INPRR3D'][1:-1,1:-1]
except KeyError:
  print 'INPRR3D not found... imposed at 0!'
  RAIN_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))

try:
  FMU_MNH=file_RSTRT.variables['FMU'][1:-1,1:-1]
  FMV_MNH=file_RSTRT.variables['FMV'][1:-1,1:-1]
  H_MNH=file_RSTRT.variables['H'][1:-1,1:-1]
  RN_MNH=file_RSTRT.variables['RN'][1:-1,1:-1]
except KeyError:
  print 'Turbulent fluxes FMU, FMV, H and LE not found... imposed at 0!'
  FMU_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))
  FMV_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))
  H_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))
  RN_MNH=np.zeros((np.shape(LON_MNH)[0],np.shape(LON_MNH)[1]))

########################################################################

print '------------------------------------------'
print ' Creating netcdf file : rstrt_MNH.nc'

fout=netCDF4.Dataset(curdir_path+'rstrt_MNH.nc','w',format='NETCDF3_64BIT')
fout.Description='Restart file for OASIS coupling'

# ----------------------------------
# Create the dimensions of the files
# ----------------------------------
fout.createDimension ('nlon', len(LON_MNH[0,:]))
fout.createDimension ('nlat', len(LAT_MNH[:,0]))

# ----------------------------------
# Create the variables of the files
# ----------------------------------
varout=fout.createVariable('MNH_FMSU','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_FMSV','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_HEAT','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_SNET','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_EVAP','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_RAIN','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH_PRES','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH__U10','d',('nlat','nlon'),fill_value=999.)
varout=fout.createVariable('MNH__V10','d',('nlat','nlon'),fill_value=999.)

# ----------------------------------
# Write out the data arrays into the file
# ----------------------------------
fout.variables['MNH_FMSU'][:,:] = FMU_MNH[:,:]
fout.variables['MNH_FMSV'][:,:] = FMV_MNH[:,:]
fout.variables['MNH_HEAT'][:,:] = H_MNH[:,:]
fout.variables['MNH_SNET'][:,:] = RN_MNH[:,:]
fout.variables['MNH_EVAP'][:,:] = EVAP_MNH[:,:]
fout.variables['MNH_RAIN'][:,:] = RAIN_MNH[:,:]
fout.variables['MNH_PRES'][:,:] = PRES_MNH[:,:]
fout.variables['MNH__U10'][:,:] = U10_MNH[:,:]
fout.variables['MNH__V10'][:,:] = V10_MNH[:,:]

# ---------------------------------------
# close the file
# ---------------------------------------
fout.close()

print ' Closing netcdf file : rstrt_MNH.nc'
print '-----------------------------------------'
#####################################################
