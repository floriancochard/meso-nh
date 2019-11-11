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
  fichier1 = addfile("KWRAI.1.SEG01.004.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_KW78")
;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = fichier1->JPHEXT
  mdims = getfilevardimsizes(fichier1,"UT") 
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2


;=================================================;
; lecture des différents champs
;=================================================;
  zs  = fichier1->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1) ; ZS
  zhat  = fichier1->ZHAT(1:kmax+1) ; ZHAT
  yhat  = fichier1->YHAT(1:jmax+jphext-1+1) ; YHAT
  xhat  = fichier1->XHAT(1:imax+jphext-1+1) ; XHAT

  wt= fichier1->WT(0,1:kmax+1,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; WT
  vt= fichier1->VT(0,1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1) ; VT
  ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1) ; UT
  tht1= fichier1->THT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; THT
  tht1@long_name="Potential Temperature"
  tht1@units="K"
  lsthm1= fichier1->LSTHM(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
  lsthm1@long_name="Large SCale Potential Temperature"
  lsthm1@units="K"

  pabst1= fichier1->PABST(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; PABST
  pabst1@long_name="Pressure"
  pabst1@units="Pa"
  rvt1= fichier1->RVT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RVT
  rvt1@long_name="Vapor mixing ratio"
  rvt1@units="kg/kg"
  rct1= fichier1->RCT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RCT
  rct1@long_name="Cloud mixing ratio"
  rct1@units="kg/kg"
  rrt1= fichier1->RRT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RRT
  rrt1@long_name="Rain mixing ratio"
  rrt1@units="kg/kg"

  inprr= fichier1->INPRR(0,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
  inprr@long_name="Instantaneous Precipitation Rain Rate"
  inprr@units="mm/h"
  acprr= fichier1->ACPRR(0,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
  acprr@long_name="Accumulated Precipitation Rain Rate"
  acprr@units="mm"

  ;YYYYDDMM1= fichier1->DTCUR__TDATE
  SS1=fichier1->DTCUR

;=================================================;
; Récupération de la date 
;=================================================;
  ;year1=YYYYDDMM1(0)
  ;day1=YYYYDDMM1(1)
  ;min1=YYYYDDMM1(2)
  time1=SS1
;=================================================;
; On mets toutes les variables sur la grille 1 
;=================================================;

  ut1 = wrf_user_unstagger(ut,"X")
  ut1@long_name="Zonal wind"
  ut1@units="m/s"
  vt1 = wrf_user_unstagger(vt,"Y")
  vt1@long_name="Meridian wind"
  vt1@units="m/s"
  wt1 = wrf_user_unstagger(wt,"Z")
  wt1@long_name="Vertical wind"
  wt1@units="m/s"

;=================================================;
; Set some other basic resources
;=================================================;
  resmap = True
  resmap@gsnFrame = False
  resmap@gsnDraw = False
  resmap@gsnMaximize = True
  resmap@gsnPaperOrientation = "portrait" 
  resmap@gsnSpreadColors= True       
  resmap@tiYAxisString =" "
  resmap@cnFillOn= True 
  resmap@cnLinesOn= False 
  resmap@lbLabelStride = 2.   ; to reduce the number of labels on labelbar
;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
  nzh=new(kmax,double)
  do k=0,kmax-1
    nzh(k)=(zhat(k)+zhat(k+1))/2.
  end do

; Create alt (altitude des niveaux modèle)
  alt=new(dimsizes(tht1),double)
  zcoef=new((/jmax,imax/),double)
  zcoef=1.-zs/nzh(kmax-1)
  do i=0,imax-1
    do j=0,jmax-1
      alt(:,j,i) = nzh*zcoef(j,i)+zs(j,i)
    end do
  end do

;=================================================;
; TRACE COUPES HORIZONTALES
;=================================================;
  gsn_define_colormap(wks,"rainbow") ; Choose colormap

  resmap@lbLabelStride = 4.   ; to reduce the number of labels on labelbar
  opts=resmap
  opts@cnLevelSelectionMode = "ManualLevels"
  opts@cnLevelSpacingF    = 0.25
  opts@cnMinLevelValF    = 0.
  opts@cnMaxLevelValF    = 6.
  plot_inprr = gsn_csm_contour(wks,inprr,opts)
  draw(plot_inprr)
  frame(wks)

  opts@cnLevelSpacingF    = 0.005
  opts@cnMinLevelValF    = 0
  opts@cnMaxLevelValF    = 0.09
  plot_acprr = gsn_csm_contour(wks,acprr,opts)
  draw(plot_acprr)
  frame(wks)

  resmap@lbLabelStride = 4.   ; to reduce the number of labels on labelbar
  res=resmap
  res@gsnScalarContour=True
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF= 10
  res@cnMinLevelValF= 98250
  res@cnMaxLevelValF= 98350
  res@vcMinDistanceF = 0.03 
  res@vcRefLengthF=0.1  
  plot=gsn_csm_vector_scalar(wks,ut1(0,:,:),vt1(0,:,:),pabst1(0,:,:),res)
  draw(plot)
  frame(wks)

; interpolattion sur les niveaux Z=3000m et Z=5000m
  rct_plane = wrf_user_intrp3d(rct1,alt,"h",3000,0.,False)
  ut_plane = wrf_user_intrp3d(ut1,alt,"h",3000,0.,False)
  vt_plane = wrf_user_intrp3d(vt1,alt,"h",3000,0.,False)
  rct5000_plane = wrf_user_intrp3d(rct1,alt,"h",5000,0.,False)
  ut5000_plane = wrf_user_intrp3d(ut1,alt,"h",5000,0.,False)
  vt5000_plane = wrf_user_intrp3d(vt1,alt,"h",5000,0.,False)

  resmap@lbLabelStride = 12.   ; to reduce the number of labels on labelbar
  res2=resmap
  res2@gsnScalarContour=True
  res2@vcMinDistanceF = 0.04
  res2@vcRefLengthF =0.1
  res2@cnLevelSpacingF = 0.00005
  res2@cnMinLevelValF = 0
  res2@tiXAxisString = "Z=3000m"
  plot=gsn_csm_vector_scalar(wks,ut_plane(:,:),vt_plane(:,:),rct_plane(:,:),res2)
  draw(plot)
  frame(wks)

  res2@tiXAxisString = "Z=5000m"
  plot=gsn_csm_vector_scalar(wks,ut5000_plane(:,:),vt5000_plane(:,:),rct5000_plane(:,:),res2)
  draw(plot)
  frame(wks)


;============================================
;   COUPE VERTICALE
;============================================
; definition de la coupe verticale
  opt_cv=False
  plane = new(2,float)
  plane =(/3,0/) ; point de depart attention i puis j
  angle=43       ; attention 0=S-N et 90=O-E

; interpolation sur la coupe verticale
  rvt_cv=wrf_user_intrp3d(rvt1,alt,"v",plane,angle,False)
  rct_cv=wrf_user_intrp3d(rct1,alt,"v",plane,angle,False)
  rrt_cv=wrf_user_intrp3d(rrt1,alt,"v",plane,angle,False)
  ut_cv=wrf_user_intrp3d(ut1,alt,"v",plane,angle,False)
  vt_cv=wrf_user_intrp3d(vt1,alt,"v",plane,angle,False)
  wt_cv=wrf_user_intrp3d(wt1,alt,"v",plane,angle,False)
  tht_cv=wrf_user_intrp3d(tht1,alt,"v",plane,angle,False)
  lsthm_cv=wrf_user_intrp3d(lsthm1,alt,"v",plane,angle,False)

  rvt_cv@long_name="Vapor mixing ratio"
  rct_cv@long_name="Cloud mixing ratio"
  rrt_cv@long_name="Rain mixing ratio"


; projection du vent horizontal dans la coupe verticale
  wind_cv=ut_cv*sin(angle*3.14/180.)+vt_cv*cos(angle*3.14/180)

  opt_cv=resmap
  opt_cv@gsnScalarContour=True
  opt_cv@vcMinDistanceF = 0.03 ; espacement des fleches
  opt_cv@vcRefLengthF=0.1      ; taille des flèches
  opt_cv@cnLevelSpacingF    = 0.0005
  opt_cv@cnMinLevelValF    = 0
  plot=gsn_csm_vector_scalar(wks,wind_cv(:,:),wt_cv(:,:),rvt_cv(:,:),opt_cv)
  draw(plot)
  frame(wks)


  res_rc=resmap
  res_rc@cnFillOn               = True ; turn on color fill
  res_rc@cnLinesOn               = False ; turn off contour lines
  res_rc@cnLevelSpacingF    = 0.00005
  res_rc@cnMinLevelValF    = 0.
  res_rc@gsnLeftString =" "
  res_rc@tiXAxisString="Cloud mixing ratio (kg/kg)"
  plot_rc=gsn_csm_contour(wks,rct_cv,res_rc)

  res_rr=resmap
  res_rr@cnFillOn               = False ; turn on color fill
  res_rr@cnLinesOn               = True ; turn off contour lines
  res_rr@cnLevelSpacingF    = 0.00025
  res_rr@cnMinLevelValF    = 0.00025
  plot_rr=gsn_csm_contour(wks,rrt_cv,res_rr)

  overlay(plot_rc,plot_rr)
  draw(plot_rc)
  frame(wks)

  res_th=resmap
  res_th@cnLevelSpacingF    = 0.2
  res_th@cnMinLevelValF    = -4
  res_th@cnMaxLevelValF    = 4
  res_th@tiMainString="THT-LSTHM (K)"
  gsn_define_colormap(wks,"BlRe") ; Choose colormap
  plot_th=gsn_csm_contour(wks,tht_cv-lsthm_cv,res_th)
  draw(plot_th)
  frame(wks)

  res_wt=resmap
  res_wt@cnLevelSpacingF    = 0.5
  res_wt@cnMinLevelValF    = -4
  res_wt@cnMaxLevelValF    = 10
  res_wt@tiMainString="WT"
  plot_wt=gsn_csm_contour(wks,wt_cv,res_wt)
  draw(plot_wt)
  frame(wks)


;=======================================
; nouvelle coupe verticale
;=======================================
  plane2 = new(2,float)
  plane2 =(/13,0/) 
  angle2=0

  rvt_cv2=wrf_user_intrp3d(rvt1,alt,"v",plane2,angle2,False)
  rct_cv2=wrf_user_intrp3d(rct1,alt,"v",plane2,angle2,False)
  rrt_cv2=wrf_user_intrp3d(rrt1,alt,"v",plane2,angle2,False)
  ut_cv2=wrf_user_intrp3d(ut1,alt,"v",plane2,angle2,False)
  vt_cv2=wrf_user_intrp3d(vt1,alt,"v",plane2,angle2,False)
  wt_cv2=wrf_user_intrp3d(wt1,alt,"v",plane2,angle2,False)
  tht_cv2=wrf_user_intrp3d(tht1,alt,"v",plane2,angle2,False)
  lsthm_cv2=wrf_user_intrp3d(lsthm1,alt,"v",plane2,angle2,False)

  gsn_define_colormap(wks,"Rainbow") 

  opt_cv=resmap
  opt_cv@gsnScalarContour=True
  opt_cv@vcMinDistanceF = 0.02 ; espacement des fleches
  opt_cv@vcRefLengthF=0.15      ; taille des flèches
  opt_cv@cnLevelSpacingF    = 0.0005
  opt_cv@cnMinLevelValF    = 0
  opt_cv@tiMainString="RVT"

  plot=gsn_csm_vector_scalar(wks,vt_cv2,wt_cv2,rvt_cv2,opt_cv)
  draw(plot)
  frame(wks)

  res_rc@tiMainString="RCT + RRT"
  plot_rc2=gsn_csm_contour(wks,rct_cv2,res_rc)
  plot_rr2=gsn_csm_contour(wks,rrt_cv2,res_rr)
  overlay(plot_rc2,plot_rr2)
  draw(plot_rc2)
  frame(wks)

  gsn_define_colormap(wks,"BlRe") ; Choose colormap
  plot_th2=gsn_csm_contour(wks,tht_cv2-lsthm_cv2,res_th)
  draw(plot_th2)
  frame(wks)
  
  plot_wt2=gsn_csm_contour(wks,wt_cv2,res_wt)
  draw(plot_wt2)
  frame(wks)


end
