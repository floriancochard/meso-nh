load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "MESONHtools.ncl"

;=============================================================;
;  Example of script to plot field on model levels from MESONH 
;  netcdf file generated with lfi2cdf
;=============================================================;

begin

;=================================================;
; Open file 
;=================================================;
; The MESONH input file.  
; ----------------------
  mnh_file="16JAN.1.12B18.001.nc"
  a = addfile(mnh_file,"r")

;==================================================;
; Open the workstation and choose colormap
; For paper quality plot do not use type ncgm or eps
; Use ps or pdf or x11 (for debugging)
;==================================================;
  type = "x11"
  wks = gsn_open_wks(type,"plt_ModelLevels")
  gsn_define_colormap(wks,"rainbow+gray")

;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"THT") ; get dimension sizes
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2

;=================================================;
; Read the variables we need
;=================================================;
  lat2d = a->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lat2d@units="degrees_north"
  lon2d = a->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lon2d@units="degrees_east"

  th  = a->THT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)   ; theta
  th@long_name="Potential Temperature"
  th@units = "K"
  th@lat2d=lat2d
  th@lon2d=lon2d

  qv  = a->RVT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)   ; Qv
  qv = qv*1000.
  qv@long_name="Water vapor mixing ratio"
  qv@units = "g/kg"
  qv@lat2d=lat2d
  qv@lon2d=lon2d

  ut   = a->UT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1)        ; u
  ut@long_name="U"
  ut@units="m/s"
  ut@lat2d=lat2d
  ut@lon2d=lon2d

  vt   = a->VT(1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1)        ; v
  vt@long_name="V"
  vt@units="m/s"
  vt@lat2d=lat2d
  vt@lon2d=lon2d

 ; Unstagger U
; --------------------
    ut!0="Z"
    ut!1="Y"
    ut!2="X"
    u = wrf_user_unstagger(ut,"X")

 ; Unstagger V
; --------------------
    vt!0="Z"
    vt!1="Y"
    vt!2="X"
    v = wrf_user_unstagger(vt,"Y")
    v@description="V"

; Calculate wind speed
; --------------------
    spd = (u*u + v*v)^(0.5)             ; speed in m/sec
    spd@long_name = "Wind Speed"
    spd@units = "m/s"

; Read projection parameters
; --------------------
    RPK  = a->RPK
    BETA = a->BETA
    LON0 = a->LON0


; Set resource for the map projection
;-----------------------------------------------;
  resmap=True

; Set map projection ressources using projection parameters
;-----------------------------------------------;
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
; Set some other basic resources
;=================================================;
  resmap@mpDataBaseVersion      = "HighRes"  ; highres map data version 
  resmap@mpFillOn               = False
  resmap@mpGeophysicalLineColor = "Black"
  resmap@mpNationalLineColor    = "Black"

  resmap@pmTickMarkDisplayMode = "Always"         ; turn on tickmarks

;=========================================================;
;  Loop over levels by step of 20
;=========================================================;
    do level =0,kmax-1,20      ; LOOP OVER LEVELS

       display_level          = level + 1
       resmap@tiMainString    = "Model Level  " + display_level

;=======================
; Resources for Theta
;=======================
       opts                                = resmap

       opts@gsnFrame              = False
       opts@gsnDraw               = False
       opts@gsnContourLineThicknessesScale = 2.0

       opts@cnLineColor                    = "Red"
       opts@cnInfoLabelOrthogonalPosF      = -0.05  
       opts@cnInfoLabelParallelPosF        = 0.97  

       contour = gsn_csm_contour_map(wks,th(level,:,:),opts)
       draw(contour)
       frame(wks)
       delete(opts)
  
;======================
; Resources for Qv
;======================
       opts = resmap

; General resources
;-------------------
       opts@gsnFrame              = False
       opts@gsnDraw               = False
       opts@gsnSpreadColors       = True
       opts@gsnSpreadColorStart   = 30
       opts@gsnSpreadColorEnd     = -2

; Contour resources
;-------------------
       opts@cnLinesOn         = False   ; no line
       opts@cnFillOn          = True    ; filled

       contour = gsn_csm_contour_map(wks,qv(level,:,:),opts)
       draw(contour)
       frame(wks)
       delete(opts)
   
;===============================================
; Overlay Wind Vectors and Speed over basic map
;===============================================

;=======================
; 1.Create Base map plot
;=======================
  resmap@gsnDraw                = False         ; don't draw yet
  resmap@gsnFrame               = False         ; don't advance frame yet

  map=gsn_map(wks,"Stereographic",resmap)

;=======================
; 2.Create Wind speed plot
;=======================
  opts_ws = True

; General resources
;-------------------
  opts_ws@gsnFrame              = False
  opts_ws@gsnDraw               = False
  opts_ws@gsnSpreadColors       = True
  opts_ws@gsnSpreadColorStart   = 30
  opts_ws@gsnSpreadColorEnd     = -2

; Tickmark resources
;---------------------
; No tickmark and label on axis for this plot
; They are already drawn by the basic map plot
; -------------------------------------------
  opts_ws@tmXBOn            = False
  opts_ws@tmYROn            = False
  opts_ws@tmXBLabelsOn      = False
  opts_ws@tmYRLabelsOn      = False

  opts_ws@lbOrientation     = "Vertical"


; Transformation resources
;---------------------------
; Field is already on the right projection
;------------------------------------------
  opts_ws@tfDoNDCOverlay = True

; Contour resources
; -------------------
  opts_ws@cnFillOn     = True
  opts_ws@cnLinesOn    = False

; PLOT
;-----
  contour = gsn_csm_contour(wks,spd(level,:,:),opts_ws)

;===========================
; 3.Create Wind vector plot
;===========================
  opts_vec = True

  opts_vec@gsnFrame   = False
  opts_vec@gsnDraw    = False
  opts_vec@gsnLeftString = "Wind Vector"
  opts_vec@gsnLeftStringParallelPosF = 0.3 
  opts_vec@gsnLeftStringOrthogonalPosF = 0.018 

; No Tickmarks/labels
;-------------
  opts_vec@tmXBOn            = False
  opts_vec@tmYROn            = False
  opts_vec@tmXBLabelsOn      = False
  opts_vec@tmYRLabelsOn      = False

  opts_vec@tfDoNDCOverlay = True

; Vector resources
;-------------------
  opts_vec@vcRefLengthF      = 0.1
  opts_vec@vcRefMagnitudeF   = 20
  opts_vec@vcMinDistanceF    = 0.05


; PLOT
;-----
   vector =  gsn_csm_vector(wks,u(level,:,:),v(level,:,:),opts_vec)

;===================
; 4.Overlay and draw
;===================
   overlay(map,contour)
   overlay(map,vector)

   draw(map)
   frame(wks)

 end do      ; END OF LEVEL LOOP
   
end
