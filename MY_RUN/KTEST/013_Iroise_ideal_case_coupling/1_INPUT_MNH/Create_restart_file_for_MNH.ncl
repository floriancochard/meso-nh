;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ================================================;
;##############################################
;#============================================#
;#      Creating restart file from prep real
;#                 case files : mesonh model
;#      Adapation from Create_restart_file_for_MNH.py
;#      Author : S. Bielli
;#      Date   :      2016
;#============================================#
;##############################################
;
begin
;
curdir_path=systemfunc("pwd")+"/"
;
file_RSTRT = curdir_path+"IROISE_5KM_201109_02_00.nc"
a=addfile(file_RSTRT,"r")

; fichier *.nc
dims = getfilevardimsizes(a,"LON")
print(dims)
imax=dims(1)-2
jmax=dims(0)-2
LON_MNH=a->LON(1:jmax,1:imax)
LAT_MNH=a->LAT(1:jmax,1:imax)
U10_MNH=a->UT(0,2,1:jmax,1:imax)
V10_MNH=a->VT(0,2,1:jmax,1:imax)
PRES_MNH=a->PABST(0,2,1:jmax,1:imax)

if(isfilevar(a,"EVAP3D")) then
  EVAP_MNH=a->EVAP3D(1:jmax,1:imax)
else
  print("EVAP3D not found... imposed at 0!")
  EVAP_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  EVAP_MNH=0.
end if

if(isfilevar(a,"INPRR3D")) then
  RAIN_MNH=a->INPRR3D(1:jmax,1:imax)
else
  print("INPRR3D not found... imposed at 0!")
  RAIN_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  RAIN_MNH=0.
end if

if(isfilevar(a,"FMU")) then
  FMU_MNH=a->FMU(1:jmax,1:imax)
  FMV_MNH=a->FMV(1:jmax,1:imax)
  H_MNH=a->H(1:jmax,1:imax)
  RN_MNH=a->RN(1:jmax,1:imax)
else
  print("Turbulent fluxes FMU, FMV, H and LE not found... imposed at 0!")
  FMU_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  FMU_MNH=0.
  FMV_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  FMV_MNH=0.
  H_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  H_MNH=0.
  RN_MNH=new(dimsizes(LON_MNH),typeof(LON_MNH))
  RN_MNH=0.
end if

;########################################################################

print("------------------------------------------")
print(" Creating netcdf file : rstrt_MNH.nc")

fout=curdir_path+"rstrt_MNH.nc"
c=addfile(fout,"c")
c@Description="Restart file for OASIS coupling"

; ----------------------------------
; Create the variables of the files
; ----------------------------------
FMU_MNH!0="nlat"
FMU_MNH!1="nlon"
FMU_MNH@_FillValue=999.

FMV_MNH!0="nlat"
FMV_MNH!1="nlon"
FMV_MNH@_FillValue=999.

H_MNH!0="nlat"
H_MNH!1="nlon"
H_MNH@_FillValue=999.

RN_MNH!0="nlat"
RN_MNH!1="nlon"
RN_MNH@_FillValue=999.

EVAP_MNH!0="nlat"
EVAP_MNH!1="nlon"
EVAP_MNH@_FillValue=999.

RAIN_MNH!0="nlat"
RAIN_MNH!1="nlon"
RAIN_MNH@_FillValue=999.

PRES_MNH!0="nlat"
PRES_MNH!1="nlon"
PRES_MNH@_FillValue=999.

U10_MNH!0="nlat"
U10_MNH!1="nlon"
U10_MNH@_FillValue=999.

V10_MNH!0="nlat"
V10_MNH!1="nlon"
V10_MNH@_FillValue=999.

; ----------------------------------
; Write out the data arrays into the file
; ----------------------------------
c->MNH_FMSU=FMU_MNH
c->MNH_FMSV=FMV_MNH
c->MNH_HEAT=H_MNH
c->MNH_SNET=RN_MNH
c->MNH_EVAP=EVAP_MNH
c->MNH_RAIN=RAIN_MNH
c->MNH_PRES=PRES_MNH
c->MNH__U10=U10_MNH
c->MNH__V10=V10_MNH

print("Closing netcdf file : rstrt_MNH.nc")
print("-----------------------------------------")
;#####################################################

end
