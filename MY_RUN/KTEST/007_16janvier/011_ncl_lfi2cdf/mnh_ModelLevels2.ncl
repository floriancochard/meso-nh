load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
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
  mnh_file="16JAN.1.12B18.001.cdf"
  a = addfile(mnh_file,"r")

;==================================================;
; Open the workstation and choose colormap
; For paper quality plot do not use type ncgm or eps
; Use ps or pdf or x11 (for debugging)
;==================================================;
  type = "x11"
  wks = gsn_open_wks(type,"plt_ModelLevels2")
  gsn_define_colormap(wks,"rainbow+gray")

;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  mdims = getfilevardimsizes(a,"THT") ; get dimension sizes
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2
  jmax=mdims(nd-2)-2
  kmax=mdims(nd-3)-2

;=================================================;
; Read the variables we need
;=================================================;
  th  = a->THT(1:kmax,1:jmax,1:imax)   ; theta
  th@long_name="Potential Temperature"
  th@units = "K"

  qv  = a->RVT(1:kmax,1:jmax,1:imax)   ; Qv
  qv = qv*1000.
  qv@long_name="Water vapor mixing ratio"
  qv@units = "g/kg"

  um   = a->UT(1:kmax,1:jmax,1:imax+1)        ; u
  um@long_name="U"
  um@units="m/s"

  vm   = a->VT(1:kmax,1:jmax+1,1:imax)        ; v
  vm@long_name="V"
  vm@units="m/s"

 ; Unstagger U
; --------------------
    um!0="Z"
    um!1="Y"
    um!2="X"
    u = wrf_user_unstagger(um,"X")

 ; Unstagger V
; --------------------
    vm!0="Z"
    vm!1="Y"
    vm!2="X"
    v = wrf_user_unstagger(vm,"Y")
    v@description="V"

; Calculate wind speed
; --------------------
    spd = (u*u + v*v)^(0.5)             ; speed in m/sec
    spd@long_name = "Wind Speed"
    spd@units = "m/s"

;=================================================;
; Use MESONHtools procedure to set map projection
; and return latitude and longitude values
;=================================================;
; Define lat and lon before calling procedure MESONH_map_c
; It is also important to define correctly the units
;----------------------------------------------------------;
  lat2d = new((/jmax,imax/),"double")
  lat2d@units="degrees_north"
  lat2d=0.
  lon2d = new((/jmax,imax/),"double")
  lon2d@units="degrees_east"
  lon2d=0.
  icornerMap = new((/2,2/),"integer")
  icornerMap(:,:)=0

; Set resource for the map projection
;-----------------------------------------------;
  resmap=True

; Call MESONH_map_c to define automatically map projection parameters
; And return lat2d, lon2d and cormap (corners) variables
;-----------------------------------------------------------------;
  MESONH_map_c (a, resmap, lat2d, lon2d, icornerMap)

;=================================================;
; Set some other basic resources
;=================================================;
  resmap@mpDataBaseVersion      = "HighRes"  ; highres map data version 
  resmap@mpFillOn               = False
  resmap@mpGeophysicalLineColor = "Black"
  resmap@mpNationalLineColor    = "Black"


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
