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
  fichier1 = addfile("REUNI.1.00A20.004.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
   type = "png"
  wks = gsn_open_wks(type,"visu_Reunion")
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
; Set map projection
;=================================================;
  lat2d=fichier1->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lon2d=fichier1->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)

; Resources for the map projection
;-----------------------------------------------;
  resmap=True
;-----------------------------------------------;
; Get global attributes to set map projection
;-----------------------------------------------;
  RPK  = fichier1->RPK ; ZS
  LON0 = fichier1->LON0
  BETA = fichier1->BETA
  if (RPK.gt.0)
    ;   Lambert  projection from north pole
    ; ---------------------------
    resmap@mpProjection          = "LambertConformal"     ; projection
    resmap@mpLambertParallel1F   =  42  ; to be adjusted !!
    resmap@mpLambertParallel2F   = resmap@mpLambertParallel1F  ;  Parallel1 = PArallel2
    resmap@mpLambertMeridianF    =  LON0      ; get value from globla attribute
  end if

  if (RPK.lt.0)
    ;   Lambert projection from south pole
    ; ---------------------------
    resmap@mpProjection          = "LambertConformal"     ; projection
    resmap@mpLambertParallel1F   =  -42  ; to be adjusted !!
    resmap@mpLambertParallel2F   = resmap@mpLambertParallel1F  ; ncl adds from grib file
    resmap@mpLambertMeridianF    =  LON0      ; get value from globla attribute
  end if

  if (RPK.eq.1)
    ;   Stereographic projection
    ; ---------------------------
    resmap@mpProjection = "Stereographic"
    resmap@mpCenterLonF = LON0
    resmap@mpCenterRotF = BETA
    resmap@mpCenterLatF = 90.
  end if

  if (RPK.eq.-1)
    ;   Stereographic projection
    ; ---------------------------
    resmap@mpProjection = "Stereographic"
    resmap@mpCenterLonF = LON0
    resmap@mpCenterRotF = BETA
    resmap@mpCenterLatF = -90.
  end if

  if (RPK.eq.0) then
    ;   Mercator projection
    ; ---------------------------
    resmap@mpProjection = "Mercator"
  end if

 print("Map projection="+resmap@mpProjection)

;====================================
; Defining the corners for projection
;====================================
  resmap@mpLimitMode            = "Corners"
  resmap@mpLeftCornerLatF       = lat2d(0,0)
  resmap@mpLeftCornerLonF       = lon2d(0,0)
  resmap@mpRightCornerLatF     = lat2d(jmax-1,imax-1)
  resmap@mpRightCornerLonF     = lon2d(jmax-1,imax-1)

;=================================================;
; lecture des différents champs
;=================================================;
  zs  = fichier1->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1) ; ZS
  zs@lat2d = lat2d
  zs@lon2d = lon2d

  wt= fichier1->WT(0,1:kmax+1,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; WT
  vt= fichier1->VT(0,1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1); VT
  ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1) ; UT
  tht= fichier1->THT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; THT
  tht@long_name="Potential Temperature"
  tht@units="K"
  tht@lat2d = lat2d
  tht@lon2d = lon2d

  lsthm= fichier1->LSTHM(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ;LSTHM
  lsthm@long_name="Large SCale Potential Temperature"
  lsthm@units="K"
  lsthm@lat2d = lat2d
  lsthm@lon2d = lon2d

  lsvm= fichier1->LSVM(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; LSVM
  lsvm@long_name="Large SCale Merdian Wind"
  lsvm@units="m/s"
  lsvm@lat2d = lat2d
  lsvm@lon2d = lon2d

  pabst= fichier1->PABST(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; PABST
  pabst@long_name="Pressure"
  pabst@units="Pa"
  pabst@lat2d = lat2d
  pabst@lon2d = lon2d


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

  ut1@lat2d = lat2d
  ut1@lon2d = lon2d
  vt1@lat2d = lat2d
  vt1@lon2d = lon2d
  wt1@lat2d = lat2d
  wt1@lon2d = lon2d
;=================================================;
; On calcule l'altitude des champs modèle
;=================================================;

  zhat= fichier1->ZHAT(1:kmax+1)

; Unstagger zhat (from grid 4 to 1)
  nzhat=new(kmax,double)
  do k=0,kmax-1
    nzhat(k)=(zhat(k)+zhat(k+1))/2.
  end do

; Create Z3D == ALT
  alt=new(dimsizes(tht),double)
  zcoef=1.-zs/nzhat(kmax-1)

  do i=0,imax-1
    do j=0,jmax-1
      alt(:,j,i) = nzhat*zcoef(j,i)+zs(j,i)
    end do
  end do

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
  resmap@lbLabelStride = 3.
  resmap@tmXBLabelStride = 3 

; If there is an error on HighRes, it means that you don't have the HighRes data
; You need to download them or the change HighRes by MediumRes
; See https://www.ncl.ucar.edu/Document/Graphics/rangs.shtml for info
; -------------------------------------
  resmap@mpDataBaseVersion    	= "MediumRes"	; choose highres map data version (must be donwloaded)

  gsn_define_colormap(wks,"rainbow") 

;=================================================;
; TRACE
;=================================================;

; module du vent + vecteur vent pour k=2

; calcul module du vent
  mod_utvt=sqrt(ut1(1,:,:)*ut1(1,:,:)+vt1(1,:,:)*vt1(1,:,:))

  mod_utvt@long_name="Module of horizontal wind"
  mod_utvt@lat2d = lat2d
  mod_utvt@lon2d = lon2d
  mod_utvt@units="m/s"

; options de tracé
  res2=resmap
  res2@cnLevelSelectionMode = "ManualLevels"
  res2@cnLevelSpacingF    =0.5
  res2@cnMinLevelValF    = 11
  res2@cnMaxLevelValF    = 26
  res2@gsnScalarContour=True
  res2@vcMinDistanceF = 0.04 
  res2@vcRefLengthF=0.1      
  res2@mpFillOn              =  False          ; turn off map fill
;  res2@mpOutlineDrawOrder    = "PostDraw"      ; draw continental outline last

  print("plot1")
  plot=gsn_csm_vector_scalar_map(wks,ut1(1,:,:),vt1(1,:,:),mod_utvt(:,:),res2)
  draw(plot)
  frame(wks)

;  temperature potentielle a 1500m
  tht_1500m = wrf_user_intrp3d(tht,alt,"h",1500,0.,False)
  res3=resmap
  res3@cnLevelSelectionMode = "ManualLevels"
  res3@cnLevelSpacingF    =0.1
  res3@cnMinLevelValF    = 302.4
  res3@cnMaxLevelValF    = 304

  print("plot2")
  plot_tht = gsn_csm_contour_map(wks,tht_1500m,res3)
  draw(plot_tht)
  frame(wks)


;  pression +vent  1500m
  pabst_1500m = wrf_user_intrp3d(pabst,alt,"h",1500,0.,False)
  ut_1500m = wrf_user_intrp3d(ut1,alt,"h",1500,0.,False)
  vt_1500m = wrf_user_intrp3d(vt1,alt,"h",1500,0.,False)

  res4=resmap
  res4@cnLevelSelectionMode = "ManualLevels"
  res4@cnLevelSpacingF=10
  res4@cnMinLevelValF= 83800
  res4@cnMaxLevelValF= 84030
  res4@gsnScalarContour=True
  res4@vcMinDistanceF = 0.04 
  res4@vcRefLengthF=0.1      
  res4@mpFillOn=  False          
  res4@mpOutlineDrawOrder= "PostDraw"     
  res4@vcMinDistanceF = 0.03 
  res4@vcRefLengthF=0.03     

  print("plot3")
  plot_pabst=gsn_csm_vector_scalar_map(wks,ut_1500m,vt_1500m,pabst_1500m,res4)
  draw(plot_pabst)
  frame(wks)

;=================================================;
;  Definition  du plan de la coupe verticale
;=================================================;
  opts=False
  plane = new(2,float)
  plane =(/33,0/)
  angle=0 
  X_plane=wrf_user_intrp2d(lat2d,plane,angle,opts)
  X_desc="latitude"

;=================================================;
; Coupe verticale de  THT
;=================================================;
; inteprolation du champ sur la coupe verticale
  tht_plane=wrf_user_intrp3d(tht,alt,"v",plane,angle,opts)
  lsthm_plane=wrf_user_intrp3d(lsthm,alt,"v",plane,angle,opts)
  vt_plane=wrf_user_intrp3d(vt1,alt,"v",plane,angle,opts)
  wt_plane=wrf_user_intrp3d(wt1,alt,"v",plane,angle,opts)
  lsvm_plane=wrf_user_intrp3d(lsvm,alt,"v",plane,angle,opts)

  ;opts_cv=resmap
  opts_cv                         = True
  opts_cv@cnLinesOn               = False
  opts_cv@cnFillOn                = True
  opts_cv@gsnSpreadColors         = True
  opts_cv@gsnFrame                = False
  opts_cv@gsnDraw                 = False
  opts_cv@lbLabelStride           = 5.

  opts_cv@tiXAxisString           = X_desc
  opts_cv@tiYAxisString           = "Height (km)"
  opts_cv@tiXAxisFontHeightF      = 0.020
  opts_cv@tiYAxisFontHeightF      = 0.020
  opts_cv@cnMissingValPerimOn     = True
  opts_cv@cnMissingValFillColor   = "gray"
  opts_cv@cnMissingValFillPattern = 17

  opts_cv@cnLevelSelectionMode = "ManualLevels"
  opts_cv@cnLevelSpacingF    =2.5
  opts_cv@cnMinLevelValF    = 300
  opts_cv@cnMaxLevelValF    = 355
  opts_cv@tiMainString="THT"
  cv_tht = gsn_csm_contour(wks,tht_plane(:,:),opts_cv)
  draw(cv_tht)
  frame(wks)

  opts_cv@cnLevelSelectionMode = "ManualLevels"
  opts_cv@cnLevelSpacingF    =0.1
  opts_cv@cnMinLevelValF    = -1
  opts_cv@cnMaxLevelValF    = 1
  opts_cv@tiMainString="THT-LSTHM"
  print("plot4")
  cv_tht_lsthm = gsn_csm_contour(wks,tht_plane(:,:)-lsthm_plane(:,:),opts_cv)
  draw(cv_tht_lsthm)
  frame(wks)

  opts_cv@cnLevelSelectionMode = "ManualLevels"
  opts_cv@cnLevelSpacingF    =0.5
  opts_cv@cnMinLevelValF    = -5
  opts_cv@cnMaxLevelValF    = 6
  opts_cv@tiMainString="VT-LSVM"
  print("plot5")
  cv_vt = gsn_csm_contour(wks,vt_plane(:,:)-lsvm_plane(:,:),opts_cv)
  draw(cv_vt)
  frame(wks)

  opts_cv@cnLevelSelectionMode = "ManualLevels"
  opts_cv@cnLevelSpacingF    =0.5
  opts_cv@cnMinLevelValF    = -6
  opts_cv@cnMaxLevelValF    = 9
  opts_cv@tiMainString="WT"
  cv_wt = gsn_csm_contour(wks,wt_plane(:,:),opts_cv)
  draw(cv_wt)
  frame(wks)

end
