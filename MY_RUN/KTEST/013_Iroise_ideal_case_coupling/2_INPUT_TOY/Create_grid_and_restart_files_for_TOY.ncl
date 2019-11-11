;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;================================================;
;###################################################
;=================================================#
;  Creating grid and restart file for toy model
;  Adaptation From Create_grid_and_restart_files_for_TOY.py
;  Author : S. Bielli
;  Date   :      2016
;=================================================#
;##################################################

begin

curdir_path=systemfunc("pwd")+"/"

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; To be defined by the user
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

;-- Limit of the grid (from etopo)
lat_domain=(/47.0, 49.5/)
lon_domain=(/-6.2, -4.0/)

;-- Type of forcing to create the restart file 
;   for the toy : CNSTE or SINUS

; CNSTE
value_CNSTE=290.0

; SINUS
value_SINUS_COEF=0.011
value_SINUS_LENGTH=1000.

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

file_topo = "topo.nc"
a=addfile(file_topo,"r")

;------ Read variables
lon_full=a->lon
lat_full=a->lat
topo_full=a->topo

ind_min_lon=ind(lon_full.gt.lon_domain(0))  
print(ind_min_lon(0))
ind_max_lon=ind(lon_full.gt.lon_domain(1))-1 
print(ind_max_lon(0))
ind_min_lat=ind(lat_full.gt.lat_domain(0)) 
print(ind_min_lat(0))
ind_max_lat=ind(lat_full.gt.lat_domain(1))-1
print(ind_max_lat(0))


lon=lon_full(ind_min_lon(0):ind_max_lon(0))
print(dimsizes(lon))
lat=lat_full(ind_min_lat(0):ind_max_lat(0))
topo=topo_full(ind_min_lat(0):ind_max_lat(0),ind_min_lon(0):ind_max_lon(0))

dims=dimsizes(topo)
nlon=dims(1) 
print("nlon="+nlon)
nlat=dims(0) 
print("nlat="+nlat)
ncorn=4           ;  print 'ncorn=', ncorn

print("---- longitude/latitude")

printVarSummary(lat)
printVarSummary(lon)
print(nlat)
print(nlon)
lat2d=conform_dims((/nlat,nlon/),lat,0)
lon2d=conform_dims((/nlat,nlon/),lon,1)
printVarSummary(lon2d)

print("---- corners longitude/latitude")
clo=new((/ncorn,nlat,nlon/),typeof(lon))
cla=new((/ncorn,nlat,nlon/),typeof(lon))
;
deltax=lon(1)-lon(0) 
print("deltax="+deltax)
clo(0,:,:)=lon2d(:,:)+deltax/2.0
clo(1,:,:)=lon2d(:,:)-deltax/2.0
clo(2,:,:)=lon2d(:,:)-deltax/2.0
clo(3,:,:)=lon2d(:,:)+deltax/2.0
;
deltay=lat(1)-lat(0) 
print("deltay="+deltay)
cla(0,:,:)=lat2d(:,:)+deltay/2.0
cla(1,:,:)=lat2d(:,:)+deltay/2.0
cla(2,:,:)=lat2d(:,:)-deltay/2.0
cla(3,:,:)=lat2d(:,:)-deltay/2.0
;
print("---- surface")
surface=new((/nlat,nlon/),typeof(lon))
surface=deltax*deltay

;
print("---- mask and var send by toy")
maskt=new((/nlat,nlon/),typeof(lon))
toyvarcnste=new((/nlat,nlon/),typeof(lon))
toyvarsinus=new((/nlat,nlon/),typeof(lon))

maskt= where(topo.gt.0.,0,1)
toyvarcnste = value_CNSTE

pi=3.141592653589793
do ind_lat=0,dimsizes(lat)-1
      toyvarsinus(ind_lat,:) = value_SINUS_COEF*sin(lat(ind_lat)*pi/180.0*value_SINUS_LENGTH)
end do

;##################################################
print("------------------------------------------")
print(" Creating netcdf file : grid_toy_model.nc")

grid_file=curdir_path+"grid_toy_model.nc"
c=addfile(grid_file,"c")
c@Description="Grid file for OASIS coupling"

; ----------------------------------
; Create the variables of the files
; ----------------------------------

lon2d!0="nlat"
lon2d!1="nlon"

lat2d!0="nlat"
lat2d!1="nlon"

clo!0="ncorner"
clo!1="nlat"
clo!2="nlon"

cla!0="ncorner"
cla!1="nlat"
cla!2="nlon"

surface!0="nlat"
surface!1="nlon"

maskt!0="nlat"
maskt!1="nlon"

;# ---------------------------------------
;# Write out the data arrays into the file
;# ---------------------------------------
c->lon=lon2d
c->lat=lat2d
c->clo=clo
c->cla=cla
c->srf=surface
c->imask=maskt

print(" Closing netcdf file : grid_toy_model.nc")
print("------------------------------------------")
;##################################################

;##################################################
print("------------------------------------------")
print(" Creating netcdf file : rstrt_TOY.nc")

rstrt_file=curdir_path+"rstrt_TOY.nc"
d=addfile(rstrt_file,"c")
d@Description="Restart file for OASIS coupling"

; ----------------------------------
; Create the variables of the files
; ----------------------------------

toyvarcnste!0="nlat"
toyvarcnste!1="nlon"

toyvarsinus!0="nlat"
toyvarsinus!1="nlon"

; ---------------------------------------
; Write out the data arrays into the file
; ---------------------------------------

d->VARCNSTE=toyvarcnste
d->VARSIN01=toyvarsinus
d->VARSIN02=toyvarsinus

print(" Closing netcdf file : rstrt_TOY.nc")
print("-----------------------------------------")
;#####################################################

end
