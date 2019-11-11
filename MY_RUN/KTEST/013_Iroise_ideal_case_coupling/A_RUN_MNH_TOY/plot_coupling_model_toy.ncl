;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;================================================;
;###################################################
;==================================================#
;      Visualization of the coupling outputs
;      Adaptation from plot_coupling_model_toy.py
;      Author : S. Bielli
;      Date   : 2016
;==================================================#
;###################################################

begin

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; To be defined by the user
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;
print("===========================================")

curdir_path=systemfunc("pwd")+"/"
print(curdir_path)

name_file_VAR1_TOY=systemfunc("ls *toyexe_01.nc")
a1=addfile(name_file_VAR1_TOY,"r")
VAR1_T=getfilevarnames(a1)
name_file_VAR1_MOD=systemfunc("ls *mesonh_01.nc")
b1=addfile(name_file_VAR1_MOD,"r")
VAR1_M=getfilevarnames(b1)
VAR1_all=stringtochar(name_file_VAR1_MOD)
VAR1=VAR1_all(5:7)
 
name_file_VAR2_TOY=systemfunc("ls *toyexe_02.nc")
a2=addfile(name_file_VAR2_TOY,"r")
VAR2_T=getfilevarnames(a2)
name_file_VAR2_MOD=systemfunc("ls *mesonh_02.nc")
b2=addfile(name_file_VAR2_MOD,"r")
VAR2_M=getfilevarnames(b2)
VAR2_all=stringtochar(name_file_VAR2_MOD)
VAR2=VAR2_all(5:7)

name_file_GRIDS="grids.nc"
g=addfile(name_file_GRIDS,"r")

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print("=============================================")
system("mkdir -p " +curdir_path+VAR1+"_"+VAR2)
print("Making directory")

print("=============================================")
print("~~~~ Plot these variables : "+ VAR1+" "+ VAR2)

;~~~~~ TOY
LAT_TOY=g->$"toyt.lat"$
LON_TOY=g->$"toyt.lon"$
dims_toy=dimsizes(LAT_TOY)
DIM_LAT_TOY=dims_toy(0)
DIM_LON_TOY=dims_toy(1)
print("DIM_LAT_TOY="+ DIM_LON_TOY+ " DIM_LAT_TOY="+DIM_LON_TOY)

;~~~~~ MOD
LAT_MOD=g->$"ssea.lat"$
LON_MOD=g->$"ssea.lon"$
dims_mod=dimsizes(LAT_MOD)
DIM_LAT_MOD=dims_mod(0)
DIM_LON_MOD=dims_mod(1)
print("DIM_LAT_MOD = "+ DIM_LON_MOD+ " DIM_LAT_MOD="+ DIM_LON_MOD)

;~~~~~ VAR1/VAR2
if (VAR1_T(0).eq."time") then
VAR1_TOY=a1->$VAR1_T(1)$
VAR1_MOD=b1->$VAR1_M(1)$
else
VAR1_TOY=a1->$VAR1_T(0)$
VAR1_MOD=b1->$VAR1_M(0)$
end if

VAR1_MOD@_FillValue = 999.
VAR1_TOY@_FillValue = 999.
VAR1_MOD = where (VAR1_MOD.gt.1000000 , VAR1_MOD@_FillValue, VAR1_MOD)
VAR1_TOY = where (VAR1_TOY.gt.1000000 , VAR1_TOY@_FillValue, VAR1_TOY)

if (VAR2_T(0).eq."time") then
VAR2_TOY=a2->$VAR2_T(1)$ 
VAR2_MOD=b2->$VAR2_M(1)$
else
VAR2_TOY=a2->$VAR2_T(0)$ 
VAR2_MOD=b2->$VAR2_M(0)$
end if

VAR1_MOD@lat2d=LAT_MOD
VAR1_MOD@lon2d=LON_MOD
VAR2_MOD@lat2d=LAT_MOD
VAR2_MOD@lon2d=LON_MOD
VAR1_TOY@lat2d=LAT_TOY
VAR1_TOY@lon2d=LON_TOY
VAR2_TOY@lat2d=LAT_TOY
VAR2_TOY@lon2d=LON_TOY

;~~~~~ TIME TOY
TIME_TOY=a1->time

;~~~~~ CONVERT VARIABLES
if (VAR2.eq."CHA") then
  print("Multiply Charnock coefficient by 1000")
  printVarSummary(VAR2_TOY)
  VAR2_TOY=VAR2_TOY*1000.0
  VAR2_MOD=VAR2_MOD*1000.0
  OPVAR2="*1000.0"
else
  OPVAR2=""
end if

OPVAR1=""

;=====================================================
;=====================================================

lev_min_VAR1 =  stringtofloat(sprintf("%5.2f",min(VAR1_TOY)))
lev_max_VAR1 =  stringtofloat(sprintf("%5.2f",max(VAR1_TOY)))
lev_int_VAR1 = stringtofloat(sprintf("%4.2f",(lev_max_VAR1-lev_min_VAR1 ) / 30.))

lev_min_VAR2 =  stringtofloat(sprintf("%5.2f",min(VAR2_MOD)))
lev_max_VAR2 =  stringtofloat(sprintf("%5.2f",max(VAR2_MOD)))
lev_int_VAR2 = stringtofloat(sprintf("%4.2f",(lev_max_VAR2-lev_min_VAR2 ) / 30.))

;-----------------------------------------------------
print("~~~~ Temporal loop")
;-----------------------------------------------------
; Open workstation and define colormap
; -------------------------------------

  type  = "png" ; open a x11 window
                ; change type to png, ps or pdf to get the plot into a file
  wks  = gsn_open_wks(type,"plt_BasicMap")   ;
  gsn_define_colormap(wks,"BlueDarkRed18")


  res                = True
  res@gsnFrame       = False
  res@gsnDraw        = False

  res@lbLabelAutoStride = True

  res@tiXAxisString         = "longitude"           ; x-axis label
  res@tiYAxisString         = "latitude"   

  res@tmXBTickSpacingF= 2
  res@tmYLTickSpacingF= 1

  res@lbOrientation   = "vertical"

  res@cnFillOn      = True
  res@cnFillMode    = "RasterFill"
  res@cnLinesOn     = False

  res@mpMaxLatF     = min((/max(LAT_MOD),max(LAT_TOY)/)) ; choose subregion
  res@mpMinLatF     = max((/min(LAT_MOD),min(LAT_TOY)/))
  res@mpMaxLonF     = min((/max(LON_TOY),max(LON_MOD)/))
  res@mpMinLonF     = max((/min(LON_TOY),min(LON_MOD)/))

  res@mpOutlineOn   = False

  res@cnLevelSelectionMode = "ManualLevels"

  res@gsnSpreadColors           = True

  resP               = True
;  gsn_define_colormap(wks,"topo_15lev") ; Choose colormap

do ind_time=0,dimsizes(TIME_TOY)-2
  plot = new(4,graphic)
  print("     ~~~~ Current time :"+ TIME_TOY(ind_time))

  resP@txString = "Cumulated time " + TIME_TOY(ind_time)+"s"
  resP@tiMainFontHeightF = 0.04

  ;----------------------
  res@tiMainString = "SND MNH " + VAR1+OPVAR1
  res@cnMinLevelValF        =  lev_min_VAR1
  res@cnMaxLevelValF        = lev_max_VAR1
  res@cnLevelSpacingF       =  lev_int_VAR1
  plot(0) =gsn_csm_contour_map(wks,VAR1_MOD(ind_time,:,:),res)

  ;----------------------
  res@tiMainString = "RCV TOY " + VAR1+OPVAR1
  plot(1) =gsn_csm_contour_map(wks,VAR1_TOY(ind_time,:,:),res)

  ;----------------------
  res@tiMainString = "SND TOY " + VAR2+OPVAR2
  res@cnMinLevelValF        =  lev_min_VAR2
  res@cnMaxLevelValF        = lev_max_VAR2
  res@cnLevelSpacingF       =  lev_int_VAR2
  plot(2) =gsn_csm_contour_map(wks,VAR2_TOY(ind_time,:,:),res)
   
  ;----------------------
  res@tiMainString = "RCV MNH " + VAR2+OPVAR2
  plot(3) =gsn_csm_contour_map(wks,VAR2_MOD(ind_time,:,:),res)

  ;------------------------

  gsn_panel(wks,plot,(/2,2/),resP)

end do
print("=============================================")

do ind_time=0,dimsizes(TIME_TOY)-2
system("mv plt_BasicMap.00000"+(ind_time+1)+".png " + curdir_path+VAR1+"_"+VAR2+"/"+VAR1+"_"+VAR2+"_MOD_TOY_T"+ind_time+".png")
end do

end
