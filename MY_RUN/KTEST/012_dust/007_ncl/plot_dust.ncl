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
  fichier1 = addfile("DUST7.1.SEG02.004.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_dust")
;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = fichier1->JPHEXT  
  mdims = getfilevardimsizes(fichier1,"UT") ; get dimension sizes
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
RPK  = fichier1->RPK 
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

 print ("Corner (0,0); Lat="+resmap@mpLeftCornerLatF+ \
                    ", Lon="+resmap@mpLeftCornerLonF)
 print ("Oppos corner; Lat="+resmap@mpRightCornerLatF+ \
                     ", Lon= "+resmap@mpRightCornerLonF)










;=================================================;
; lecture des différents champs
;=================================================;
wt= fichier1->WT(0,1:kmax+1,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; WT
vt= fichier1->VT(0,1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1) ; VT
ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1) ; UT

tht= fichier1->THT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; THT
tht@long_name="Potential Temperature"
tht@units="K"
tht@lat2d = lat2d
tht@lon2d = lon2d



dstm03t= fichier1->DSTM03T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm03t@long_name="DSTM03T"
dstm03t@units="ppp"
dstm03t@lat2d = lat2d
dstm03t@lon2d = lon2d

dstm33t= fichier1->DSTM33T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm33t@long_name="DSTM33T"
dstm33t@units="ppp"
dstm33t@lat2d = lat2d
dstm33t@lon2d = lon2d

dstm02t= fichier1->DSTM02T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm02t@long_name="DSTM02T"
dstm02t@units="ppp"
dstm02t@lat2d = lat2d
dstm02t@lon2d = lon2d

dstm32t= fichier1->DSTM32T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm32t@long_name="DSTM32T"
dstm32t@units="ppp"
dstm32t@lat2d = lat2d
dstm32t@lon2d = lon2d

dstm01t= fichier1->DSTM01T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm01t@long_name="DSTM01T"
dstm01t@units="ppp"
dstm01t@lat2d = lat2d
dstm01t@lon2d = lon2d

dstm31t= fichier1->DSTM31T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
dstm31t@long_name="DSTM31T"
dstm31t@units="ppp"
dstm31t@lat2d = lat2d
dstm31t@lon2d = lon2d

flx_dst001= fichier1->F_DST001P1(jphext:jmax+jphext-1,jphext:imax+jphext-1)
flx_dst001@long_name="F_DST001P1"
flx_dst001@units="kg/m2/s"
flx_dst001@lat2d = lat2d
flx_dst001@lon2d = lon2d


flx_dst002= fichier1->F_DST002P1(jphext:jmax+jphext-1,jphext:imax+jphext-1)
flx_dst002@long_name="F_DST002P1"
flx_dst002@units="kg/m2/s"
flx_dst002@lat2d = lat2d
flx_dst002@lon2d = lon2d

flx_dst003= fichier1->F_DST003P1(jphext:jmax+jphext-1,jphext:imax+jphext-1)
flx_dst003@long_name="F_DST003P1"
flx_dst003@units="kg/m2/s"
flx_dst003@lat2d = lat2d
flx_dst003@lon2d = lon2d

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
; Set some other basic resources
;=================================================;
  resmap = True
  resmap@gsnFrame = False
  resmap@gsnDraw = False
  resmap@gsnMaximize = True
  resmap@gsnPaperOrientation = "portrait" 
  resmap@gsnSpreadColors    	= True       	; use full range of colormap
  resmap@tiYAxisString =" "
  resmap@cnFillOn               = True ; turn on color fill
  resmap@cnLinesOn               = False ; turn off contour lines
  resmap@mpFillOn               = False ; couleur grise des continents
  resmap@mpGeophysicalLineThicknessF = 2.0
  resmap@mpGeophysicalLineColor = "Black"
  resmap@mpNationalLineColor    = "Black"

;=================================================;
; module + vecteurs vent au premier niveau
;=================================================;
 
gsn_define_colormap(wks,"rainbow") ; Choose colormap

; options de tracé
res2=resmap
res2@cnLevelSelectionMode = "AutomaticLevels"
res2@gsnScalarContour=True
res2@vcMinDistanceF = 0.01 
res2@vcRefLengthF=0.04      
res2@cnLevelSpacingF    = 0.00005
res2@cnMinLevelValF    = 0

mutvt=(ut1*ut1+vt1*vt1)^(0.5)
mutvt@lat2d=lat2d
mutvt@lon2d=lon2d

plot=gsn_csm_vector_scalar_map(wks,ut1(0,:,:),vt1(0,:,:),mutvt(0,:,:),res2)
draw(plot)
frame(wks)

plot=gsn_csm_contour_map(wks,tht(0,:,:),res2)
draw(plot)
frame(wks)


res2@gsnSpreadColorStart  = 0     ; force la première couleur en blanc(= couleur 0 de la palette)
plot=gsn_csm_contour_map(wks,dstm03t(0,:,:),res2)
draw(plot)
frame(wks)
plot=gsn_csm_contour_map(wks,dstm33t(0,:,:),res2)
draw(plot)
frame(wks)

plot=gsn_csm_contour_map(wks,dstm02t(0,:,:),res2)
draw(plot)
frame(wks)
plot=gsn_csm_contour_map(wks,dstm32t(0,:,:),res2)
draw(plot)
frame(wks)

plot=gsn_csm_contour_map(wks,dstm01t(0,:,:),res2)
draw(plot)
frame(wks)
plot=gsn_csm_contour_map(wks,dstm31t(0,:,:),res2)
draw(plot)
frame(wks)


flx_dst001@_FillValue = 999.
plot=gsn_csm_contour_map(wks,flx_dst001(:,:),res2)
draw(plot)
frame(wks)
flx_dst002@_FillValue = 999.
plot=gsn_csm_contour_map(wks,flx_dst002(:,:),res2)
draw(plot)
frame(wks)
flx_dst003@_FillValue = 999.
plot=gsn_csm_contour_map(wks,flx_dst003(:,:),res2)
draw(plot)
frame(wks)


end
