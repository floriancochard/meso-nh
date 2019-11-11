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
  fichier1 = addfile("16JAN.2.12B18.001.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_16j_2")
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
zs  = fichier1->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1) ; ZS
zs@long_name="Orography"
zs@units="m"
zs@lat2d = lat2d
zs@lon2d = lon2d


wt= fichier1->WT(0,1:kmax+1,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; WT
vt= fichier1->VT(0,1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1) ; VT
ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1) ; UT

tht1= fichier1->THT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; THT
tht1@long_name="Potential Temperature"
tht1@units="K"
tht1@lat2d = lat2d
tht1@lon2d = lon2d



pabst1= fichier1->PABST(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; PABST
pabst1@long_name="Pressure"
pabst1@units="Pa"
pabst1@lat2d = lat2d
pabst1@lon2d = lon2d

rvt1= fichier1->RVT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RVT
rvt1=rvt1*1000
rvt1@long_name="Vapor mixing ratio"
rvt1@units="kg/kg"
rvt1@lat2d = lat2d
rvt1@lon2d = lon2d

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
    alt=new(dimsizes(tht1),double)
    zcoef=1.-zs/nzhat(kmax-1)

    do i=0,imax-1
    do j=0,jmax-1
       alt(:,j,i) = nzhat*zcoef(j,i)+zs(j,i)
    end do
    end do

alt@lat2d = lat2d
alt@lon2d = lon2d

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

  resmap@tmXBLabelStride = 2   ; to reduce the number of labels on xaxis
  resmap@lbLabelStride = 2.   ; to reduce the number of labels on labelbar


;=================================================;
; THT + vecteurs vent a 850 hPa
;=================================================;
 
gsn_define_colormap(wks,"rainbow") ; Choose colormap

; interpolation des champs a 850 hPa
tht_plane = wrf_user_intrp3d(tht1,pabst1,"h",850,0.,False)
ut_plane = wrf_user_intrp3d(ut1,pabst1,"h",850,0.,False)
vt_plane = wrf_user_intrp3d(vt1,pabst1,"h",850,0.,False)

; options de tracé
res2=resmap
res2@cnLevelSelectionMode = "AutomaticLevels"
res2@gsnScalarContour=True
res2@vcMinDistanceF = 0.04 
res2@vcRefLengthF=0.1      
res2@tiXAxisString = "850 hPa"


plot=gsn_csm_vector_scalar_map(wks,ut_plane(:,:),vt_plane(:,:),tht_plane(:,:),res2)
draw(plot)
frame(wks)

;=================================================;
; RVT + vecteurs vent a 700 hPa
;=================================================;

rvt_plane = wrf_user_intrp3d(rvt1,pabst1,"h",700,0.,False)

res2@tiXAxisString = "700 hPa"

plot=gsn_csm_vector_scalar_map(wks,ut_plane(:,:),vt_plane(:,:),rvt_plane(:,:),res2)
draw(plot)
frame(wks)

;=================================================;
; PABST + vecteurs vent a 9000 m
;=================================================;

pabst9000_plane = wrf_user_intrp3d(pabst1,alt,"h",9000,0.,False)
ut9000_plane = wrf_user_intrp3d(ut1,alt,"h",9000,0.,False)
vt9000_plane = wrf_user_intrp3d(vt1,alt,"h",9000,0.,False)

res2@tiXAxisString = "Z=9000m"

plot=gsn_csm_vector_scalar_map(wks,ut9000_plane,vt9000_plane,pabst9000_plane,res2)
draw(plot)
frame(wks)

end
