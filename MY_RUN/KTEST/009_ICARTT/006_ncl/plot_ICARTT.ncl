;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
; ================================================;
begin
;=================================================;
; open file and read in data
;=================================================;
  a = addfile("ICART.1.SEG01.001dg.nc", "r")
  a2 = addfile("ICART.1.SEG01.002dg.nc", "r")

;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"UT") 
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2

;-------------------------------------------------;
; Read data.
;-------------------------------------------------;
  lat2d = a->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lat2d@units="degrees_north"
  lon2d = a->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lon2d@units="degrees_east"

zs  = a->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1) ; ZS
zs@long_name="Orography"
zs@units="m"
zs@lat2d = lat2d
zs@lon2d = lon2d

printMinMax(zs,0)

  rc_t1 = a->MRC(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  rc_t1@long_name="Cloud mixing ratio"
  rc_t1@units="g/kg"
  rc_t1@lat2d=lat2d
  rc_t1@lon2d=lon2d
printMinMax(rc_t1,0)

;
  o3_t1 = a->O3T(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3_t1@long_name="Ozone"
  o3_t1@units="ppbv"
  o3_t1@lat2d=lat2d
  o3_t1@lon2d=lon2d

;
  co_t1 = a->COT(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  co_t1@long_name="carbon monoxide"
  co_t1@units="ppbv"
  co_t1@lat2d=lat2d
  co_t1@lon2d=lon2d

;
  o3p_t1=a->O3_PROD(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3p_t1 = o3p_t1*1e9*3600
  o3p_t1@long_name="ozone production"
  o3p_t1@units="ppbv/h"  
  o3p_t1@lat2d=lat2d
  o3p_t1@lon2d=lon2d  

;
  o3l_t1=a->O3_LOSS(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3l_t1 = o3l_t1*1e9*3600
  o3l_t1@long_name="ozone loss"
  o3l_t1@units="ppbv/h"  
  o3l_t1@lat2d=lat2d
  o3l_t1@lon2d=lon2d  

;
  cop_t1=a->CO_PROD(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  cop_t1 = cop_t1*1e9*3600
  cop_t1@long_name="carbon monoxide production"
  cop_t1@units="ppbv/h"  
  cop_t1@lat2d=lat2d
  cop_t1@lon2d=lon2d  

;   
  col_t1=a->CO_LOSS(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  col_t1 = col_t1*1e9*3600
  col_t1@long_name="carbon monoxide loss"
  col_t1@units="ppbv/h"
  col_t1@lat2d=lat2d
  col_t1@lon2d=lon2d


;
;
  rc_t2 = a2->MRC(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  rc_t2@long_name="Cloud mixing ratio"
  rc_t2@units="g/kg"
  rc_t2@lat2d=lat2d
  rc_t2@lon2d=lon2d

;
  o3_t2 = a2->O3T(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3_t2@long_name="Ozone"
  o3_t2@units="ppbv"
  o3_t2@lat2d=lat2d
  o3_t2@lon2d=lon2d

;
  co_t2 = a2->COT(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  co_t2@long_name="carbon monoxide"
  co_t2@units="ppbv"
  co_t2@lat2d=lat2d
  co_t2@lon2d=lon2d

;
  o3p_t2=a2->O3_PROD(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3p_t2=o3p_t2*1e9*3600
  o3p_t2@long_name="ozone production"
  o3p_t2@units="pptv/h"  
  o3p_t2@lat2d=lat2d
  o3p_t2@lon2d=lon2d  

;
  o3l_t2=a2->O3_LOSS(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3l_t2=o3l_t2*1e9*3600
  o3l_t2@long_name="ozone loss"
  o3l_t2@units="pptv/h"
  o3l_t2@lat2d=lat2d
  o3l_t2@lon2d=lon2d

;
  cop_t2=a2->CO_PROD(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  cop_t2=cop_t2*1e9*3600
  cop_t2@long_name="carbon monoxide production"
  cop_t2@units="pptv/h"
  cop_t2@lat2d=lat2d
  cop_t2@lon2d=lon2d

;   
  col_t2=a2->CO_LOSS(:,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  col_t2=col_t2*1e9*3600
  col_t2@long_name="carbon monoxide loss"
  col_t2@units="ppbv/h"
  col_t2@lat2d=lat2d
  col_t2@lon2d=lon2d  

  

;-----------------------------------------------;
;=================================================;
; On calcule l'altitude des champs modèle
;=================================================;

zhat= a2->ZHAT(1:kmax+1)

; Unstagger zhat (from grid 4 to 1)
    nzhat=new(kmax,double)
    do k=0,kmax-1
     nzhat(k)=(zhat(k)+zhat(k+1))/2.
    end do

; Create Z3D == ALT
    alt=new(dimsizes(o3_t2),double)
    zcoef=1.-zs/nzhat(kmax-1)

    do i=0,imax-1
    do j=0,jmax-1
       alt(0,:,j,i) = nzhat*zcoef(j,i)+zs(j,i)
    end do
    end do

alt@lat2d = lat2d
alt@lon2d = lon2d




;-----------------------------------------------;
; Set map projection ressources using projection parameters
;-----------------------------------------------;
; Read projection parameters
; --------------------
    RPK  = a2->RPK
    BETA = a2->BETA
    LON0 = a2->LON0

  resmap=True
  if (RPK.gt.0)
; ---------------------------
  ;   Lambert  projection from north pole
; ---------------------------
   resmap@mpProjection          = "LambertConformal"     ; projection
   pole                         = 1    ; projection for north hemisphere
   resmap@mpLambertParallel1F   = pole*asin(RPK)*180/3.14 ; with pole=1 for north hemisphere and -1 for south hemisphere
   resmap@mpLambertParallel2F   = resmap@mpLambertParallel1F  ; ncl adds from grib file
   resmap@mpLambertMeridianF    = LON0      ; ncl adds from grib file
  end if

  if (RPK.lt.0)
; ---------------------------
  ;   Lambert projection from south pole
; ---------------------------
   resmap@mpProjection          = "LambertConformal"     ; projection
   pole                         = -1                     ; projection for south hemisphere
   resmap@mpLambertParallel1F   = pole*asin(RPK)*180/3.14 ; with pole=1 for north hemisphere and -1 for south hemisphere
   resmap@mpLambertParallel2F   = resmap@mpLambertParallel1F  ; ncl adds from grib file
   resmap@mpLambertMeridianF    =  LON0      ; ncl adds from grib file
  end if

  if (RPK.eq.1)
; ---------------------------
  ;   Stereographic projection north
; ---------------------------
    resmap@mpProjection = "Stereographic"
    resmap@mpCenterLonF           = LON0
    resmap@mpCenterRotF           = BETA
    resmap@mpCenterLatF           = 90
  end if

  if (RPK.eq.-1)
; ---------------------------
  ;   Stereographic projection south
; ---------------------------
    resmap@mpProjection = "Stereographic"
    resmap@mpCenterLonF           = LON0
    resmap@mpCenterRotF           = BETA
    resmap@mpCenterLatF           = -90
  end if

  if (RPK.eq.0) then
; ---------------------------
  ;   Mercator projection
; ---------------------------
    resmap@mpProjection = "Mercator"
  end if

 print("Map projection="+resmap@mpProjection)

; Defining the corners for projection
; --------------------------------
  resmap@mpLimitMode            = "Corners"
  resmap@mpLeftCornerLatF       = lat2d(0,0)
  resmap@mpLeftCornerLonF       = lon2d(0,0)
  resmap@mpRightCornerLatF     = lat2d(jmax-1,imax-1)
  resmap@mpRightCornerLonF     = lon2d(jmax-1,imax-1)
  
;=================================================;
; PLOT
;=================================================;
; interpolation des champs a 1250 m
rc_t1_plane = wrf_user_intrp3d(rc_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
printMinMax(rc_t1_plane,0)
printMinMax(alt,0)

rc_t2_plane = wrf_user_intrp3d(rc_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
co_t1_plane = wrf_user_intrp3d(co_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
co_t2_plane = wrf_user_intrp3d(co_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
cop_t1_plane= wrf_user_intrp3d(cop_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
cop_t2_plane= wrf_user_intrp3d(cop_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
col_t1_plane= wrf_user_intrp3d(col_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
col_t2_plane= wrf_user_intrp3d(col_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3_t1_plane = wrf_user_intrp3d(o3_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3_t2_plane = wrf_user_intrp3d(o3_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3l_t1_plane= wrf_user_intrp3d(o3l_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3l_t2_plane= wrf_user_intrp3d(o3l_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3p_t1_plane= wrf_user_intrp3d(o3p_t1(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)
o3p_t2_plane= wrf_user_intrp3d(o3p_t2(0,:,:,:),alt(0,:,:,:),"h",1250,0.,False)


  figname ="zsection_1250"
  wks  = gsn_open_wks("png",figname)   ; open a ncgm file
  gsn_define_colormap(wks,"WhBlGrYeRe") ; Choose colormap

  res                 = resmap          
  res@gsnDraw   		= False         ; don't draw yet
  res@gsnFrame  		= False         ; don't advance frame yet

; X-axis title (tiY)                              
  res@tiXAxisFontHeightF = 0.018  ; font height
  res@tiXAxisFont        = 21     ; font index
  res@tiXAxisString      = "longitude"  ; string to use as the X-Axis title

; Y-axis title (tiY)
  res@tiYAxisFontHeightF = 0.018  ; font height
  res@tiYAxisFont        = 21     ; font index
  res@tiYAxisString      = "latitude"  ; string to use as the Y-Axis title

; BW
  res@cnLinesOn            = False           
  res@cnFillOn             = True        
  res@gsnSpreadColors      = True 
;
; label bar (lb)
;  res@lbAutoManage       = False  
;  res@lbBottomMarginF    = 0.4        ; offset
;  res@lbOrientation      = "Vertical"

; Map ressources 
;  res@mpDataBaseVersion    	= "HighRes"	; choose highres map data version (must be donwloaded)
;  res@mpDataBaseVersion    	= "MediumRes"	; choose highres map data version (must be donwloaded)
  res@mpGridAndLimbOn   	= True             ; turn on lat/lon lines
  res@mpGridLatSpacingF 	= 10.              ; spacing for lat lines
  res@mpGridLonSpacingF 	= 10.              ; spacing for lon lines

  res@mpGeophysicalLineColor = "Black"  ; default value in lowres
  res@mpNationalLineColor    = "Black"  ; idem
  res@mpUSStateLineColor     = "Black"  ; idem
  res@mpGridLineColor        = "Black"
  res@mpLimbLineColor        = "Black"
  res@mpPerimLineColor       = "Black"

  
  res@gsnCenterString="heure=19"

; plot cloud mixing ratio
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels =     (/0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06/)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96,101/) ; color of a contour
  plot_rc = gsn_csm_contour_map(wks,rc_t1_plane(:,:),res)
  draw(plot_rc)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot ozone
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/15.,20., 25., 35., 40., 45., 50., 55., 60., 65./)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o3 = gsn_csm_contour_map(wks,o3_t1_plane(:,:),res)
  draw(plot_o3)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot co
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/110.,112.5,115.,117.5,120.,122.5,125.,127.5,130.,132.5,135./)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96,101/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_co = gsn_csm_contour_map(wks,co_t1_plane(:,:),res)
  draw(plot_co)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot ozone production
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/75,80,85,90,95,100,105,110,115,120,125/)
  res@cnFillColors = (/75,80,85,90,95,100,105,110,115,120,125,130/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o3p = gsn_csm_contour_map(wks,o3p_t1_plane(:,:),res)
  draw(plot_o3p)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot ozone loss
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/75,80,85,90,95,100,105,110,115,120,125/)
  res@cnFillColors = (/75,80,85,90,95,100,105,110,115,120,125,130/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o3l = gsn_csm_contour_map(wks,o3l_t1_plane(:,:),res)
  draw(plot_o3l)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot carbon monoxide production
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/)
  res@cnFillColors := (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_cop = gsn_csm_contour_map(wks,cop_t1_plane(:,:),res)
  draw(plot_cop)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot carbon monoxide loss
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/)
  res@cnFillColors := (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_col = gsn_csm_contour_map(wks,col_t1_plane(:,:),res)
  draw(plot_col)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

  res@gsnCenterString="heure=20"

; plot cloud mixing ratio
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels =     (/0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06/)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96,101/) ; color of a contour
  plot_rc1 = gsn_csm_contour_map(wks,rc_t2_plane(:,:),res)
  draw(plot_rc1)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot ozone
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/15.,20., 25., 35., 40., 45., 50., 55., 60., 65./)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o31 = gsn_csm_contour_map(wks,o3_t2_plane(:,:),res)
  draw(plot_o31)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot co
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/110.,112.5,115.,117.5,120.,122.5,125.,127.5,130.,132.5,135./)
  res@cnFillColors = (/2,6,12,40,45,51,62,72,80,89,96,101/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_co1 = gsn_csm_contour_map(wks,co_t2_plane(:,:),res)
  draw(plot_co1)
  frame(wks)

; plot ozone production
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/75,80,85,90,95,100,105,110,115,120,125/)
  res@cnFillColors = (/75,80,85,90,95,100,105,110,115,120,125,130/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o3p1 = gsn_csm_contour_map(wks,o3p_t2_plane(:,:),res)
  draw(plot_o3p1)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot ozone loss
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/75,80,85,90,95,100,105,110,115,120,125/)
  res@cnFillColors = (/75,80,85,90,95,100,105,110,115,120,125,130/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_o3l1 = gsn_csm_contour_map(wks,o3l_t2_plane(:,:),res)
  draw(plot_o3l1)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot carbon monoxide production
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/)
  res@cnFillColors := (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_cop1 = gsn_csm_contour_map(wks,cop_t2_plane(:,:),res)
  draw(plot_cop1)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

; plot carbon monoxide loss
  res@cnLevelSelectionMode = "ExplicitLevels" ; method for selecting the contour
  res@cnLevels     = (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/)
  res@cnFillColors := (/0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6/) ; color of a contour
;  res@cnLevelSelectionMode = "AutomaticLevels" 
  plot_col1 = gsn_csm_contour_map(wks,col_t2_plane(:,:),res)
  draw(plot_col1)
  frame(wks)
  delete(res@cnLevels)
  delete(res@cnFillColors)

;;;;;;;;;;;;;;;;;;;;;;;;

end



