##!/usr/bin/python
# -*- coding: utf-8 -*-
#
####################################################
#==================================================#
#      Visualization of the coupling outputs
#      Author : J. Pianezze
#      Date   : 2015
#==================================================#
####################################################

import netCDF4
import numpy as np
import scipy
import matplotlib.pyplot as plt
from pylab import *
import os
import glob

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# To be defined by the user
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
print '=============================================='
curdir_path=os.path.abspath(os.curdir)+'/'
print curdir_path
#
name_file_VAR1_TOY=glob.glob('*toyexe_01.nc')[0]
name_file_VAR1_MOD=glob.glob('*mesonh_01.nc')[0]
VAR1=name_file_VAR1_MOD[5:8]
#
name_file_VAR2_TOY=glob.glob('*toyexe_02.nc')[0]
name_file_VAR2_MOD=glob.glob('*mesonh_02.nc')[0]
VAR2=name_file_VAR2_MOD[5:8]
#
name_file_GRIDS='grids.nc'
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print '=============================================='
try :
  os.mkdir(curdir_path+VAR1+'_'+VAR2+'/')
except OSError:
  print 'Directory already created'
else:
  print 'Making directory'

print '=============================================='
print '~~~~ Plot these variables :', VAR1, VAR2

file_VAR1_TOY = netCDF4.Dataset(name_file_VAR1_TOY)
file_VAR1_MOD = netCDF4.Dataset(name_file_VAR1_MOD)

file_VAR2_TOY = netCDF4.Dataset(name_file_VAR2_TOY)
file_VAR2_MOD = netCDF4.Dataset(name_file_VAR2_MOD)
 
file_GRIDS = netCDF4.Dataset(name_file_GRIDS)

#~~~~~ TOY
LAT_TOY=file_GRIDS.variables['toyt.lat']
LON_TOY=file_GRIDS.variables['toyt.lon']
DIM_LAT_TOY=np.shape(LAT_TOY)[0]
DIM_LON_TOY=np.shape(LON_TOY)[1]
print 'DIM_LAT_TOY, DIM_LON_TOY', DIM_LAT_TOY, DIM_LON_TOY

#~~~~~ MOD
LAT_MOD=file_GRIDS.variables['ssea.lat']
LON_MOD=file_GRIDS.variables['ssea.lon']
DIM_LAT_MOD=np.shape(LAT_MOD)[0]
DIM_LON_MOD=np.shape(LON_MOD)[1]
print 'DIM_LAT_MOD, DIM_LON_MOD', DIM_LAT_MOD, DIM_LON_MOD

#~~~~~ VAR1/VAR2
VAR1_TOY=file_VAR1_TOY.variables[name_file_VAR1_TOY[0:8]][:,:,:]
VAR1_MOD=file_VAR1_MOD.variables[name_file_VAR1_MOD[0:8]][:,:,:]
MASK_VAR1_MOD = (VAR1_MOD[:,:,:] == 1E20)
VAR1_MOD = np.ma.MaskedArray(VAR1_MOD, mask=MASK_VAR1_MOD)

VAR2_TOY=file_VAR2_TOY.variables[name_file_VAR2_TOY[0:8]][:,:,:] 
VAR2_MOD=file_VAR2_MOD.variables[name_file_VAR2_MOD[0:8]][:,:,:]

OPVAR1=''
OPVAR2=''

#~~~~~ TIME TOY
TIME_TOY=file_VAR1_TOY.variables['time'][:]

#~~~~~ CONVERT VARIABLES
if VAR2=='CHA':
  print 'Multiply Charnock coefficient by 1000'
  VAR2_TOY[:,:,:]=VAR2_TOY[:,:,:]*1000.0
  VAR2_MOD[:,:,:]=VAR2_MOD[:,:,:]*1000.0
  OPVAR2='*1000.0'

#=====================================================
#=====================================================

#-----------------------------------------------------
print '~~~~ Temporal loop'
#-----------------------------------------------------
for ind_time in xrange(np.size(TIME_TOY)-1):
  print '     ~~~~ Current time :', TIME_TOY[ind_time]

  fig = plt.figure()
  fig.suptitle('Cumulated time : '+str(TIME_TOY[ind_time])+'s', fontsize=18)

  #----------------------
  ax = fig.add_subplot(221)
  plt.title('SND MNH '+VAR1+OPVAR1,fontsize=18)
  CS=plt.pcolormesh(LON_MOD[:,:],LAT_MOD[:,:],VAR1_MOD[ind_time,:,:],cmap=plt.cm.RdBu_r,\
  	                vmin=np.min(VAR1_TOY), vmax=np.max(VAR1_TOY))
  cbar = plt.colorbar(CS,orientation='vertical',format='%.3f')
  plt.ylabel(r'latitude [-]',fontsize=18)
  plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
  xlim(( max(np.min(LON_MOD[:,:]),np.min(LON_TOY[:,:])), min(np.max(LON_MOD[:,:]),np.max(LON_TOY[:,:])) ))
  ylim(( max(np.min(LAT_MOD[:,:]),np.min(LAT_TOY[:,:])), min(np.max(LAT_MOD[:,:]),np.max(LAT_TOY[:,:])) ))

  #----------------------
  ax = fig.add_subplot(222)
  plt.title('RCV TOY '+VAR1+OPVAR1,fontsize=18)
  CS=plt.pcolormesh(LON_TOY[:,:],LAT_TOY[:,:],VAR1_TOY[ind_time,:,:],cmap=plt.cm.RdBu_r,\
  	                vmin=np.min(VAR1_TOY), vmax=np.max(VAR1_TOY))  
  cbar = plt.colorbar(CS,orientation='vertical',format='%.3f')
  plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
  plt.tick_params(axis='y',which='both',bottom='on',top='on',labelleft='off')
  xlim(( max(np.min(LON_MOD[:,:]),np.min(LON_TOY[:,:])), min(np.max(LON_MOD[:,:]),np.max(LON_TOY[:,:])) ))
  ylim(( max(np.min(LAT_MOD[:,:]),np.min(LAT_TOY[:,:])), min(np.max(LAT_MOD[:,:]),np.max(LAT_TOY[:,:])) ))

  #----------------------
  ax = fig.add_subplot(223)
  plt.title('SND TOY '+VAR2+OPVAR2,fontsize=18)
  CS=plt.pcolormesh(LON_TOY[:,:],LAT_TOY[:,:],VAR2_TOY[ind_time,:,:],cmap=plt.cm.RdBu_r,\
  	                vmin=np.min(VAR2_MOD), vmax=np.max(VAR2_MOD))  
  cbar = plt.colorbar(CS,orientation='vertical',format='%.3f')
  plt.ylabel(r'latitude [-]',fontsize=18)
  plt.xlabel(r'longitude [-]',fontsize=18)
  xlim(( max(np.min(LON_MOD[:,:]),np.min(LON_TOY[:,:])), min(np.max(LON_MOD[:,:]),np.max(LON_TOY[:,:])) ))
  ylim(( max(np.min(LAT_MOD[:,:]),np.min(LAT_TOY[:,:])), min(np.max(LAT_MOD[:,:]),np.max(LAT_TOY[:,:])) ))
   
  #----------------------
  ax = fig.add_subplot(224)
  plt.title('RCV MNH '+VAR2+OPVAR2,fontsize=18)
  CS=plt.pcolormesh(LON_MOD[:,:],LAT_MOD[:,:],VAR2_MOD[ind_time,:,:],cmap=plt.cm.RdBu_r,\
  	                vmin=np.min(VAR2_MOD), vmax=np.max(VAR2_MOD))  
  cbar = plt.colorbar(CS,orientation='vertical',format='%.3f')
  plt.xlabel(r'longitude [-]',fontsize=18)
  plt.tick_params(axis='y',which='both',bottom='on',top='on',labelleft='off')
  xlim(( max(np.min(LON_MOD[:,:]),np.min(LON_TOY[:,:])), min(np.max(LON_MOD[:,:]),np.max(LON_TOY[:,:])) ))
  ylim(( max(np.min(LAT_MOD[:,:]),np.min(LAT_TOY[:,:])), min(np.max(LAT_MOD[:,:]),np.max(LAT_TOY[:,:])) ))

  #------------------------
  plt.savefig(curdir_path+VAR1+"_"+VAR2+"/"+VAR1+"_"+VAR2+"_MOD_TOY_T"+str(ind_time)+".pdf")
  plt.savefig(curdir_path+VAR1+"_"+VAR2+"/"+VAR1+"_"+VAR2+"_MOD_TOY_T"+str(ind_time)+".png")
  plt.show()
  plt.close()
print '=============================================='
