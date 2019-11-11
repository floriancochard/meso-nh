;================================================;
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
; ================================================;
;================================================;
  begin
;=================================================;
; open file and read in data
;=================================================;
  fichier1 = addfile("REL3D.1.EXP01.002.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_3Drelief")
;=================================================;
; lecture des différents champs
;=================================================;
  jphext  = fichier1->JPHEXT
  mdims = getfilevardimsizes(fichier1,"UT") 
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2



  ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; UT
  ut@long_name="Zonal wind"
  ut@units="m/s"



;=================================================;
; Unstagger U
;=================================================;
  ut!0="Z"
  ut!1="Y"
  ut!2="X"
  u = wrf_user_unstagger(ut,"X")

;=================================================;
; Set some other basic resources
;=================================================;
  resmap = True
  resmap@gsnFrame = False
  resmap@gsnDraw = False
  resmap@gsnMaximize = True
  resmap@gsnPaperOrientation = "portrait" 
  resmap@gsnSpreadColors= True       	
  resmap@tiMainString ="UT level 2 "
  resmap@cnFillOn= True
  resmap@cnLinesOn= False 

  resmap@lbLabelStride = 4. ; to reduce number of labels on labelbar


;=================================================;
; TRACE
;=================================================;
  gsn_define_colormap(wks,"rainbow")

  resmap@cnLevelSelectionMode = "ManualLevels"
  resmap@cnLevelSpacingF    = 0.25
  resmap@cnMinLevelValF    = 8
  resmap@cnMaxLevelValF    = 14
  plot_ut = gsn_csm_contour(wks,u(0,:,:),resmap)
  draw(plot_ut)
  frame(wks)

end
