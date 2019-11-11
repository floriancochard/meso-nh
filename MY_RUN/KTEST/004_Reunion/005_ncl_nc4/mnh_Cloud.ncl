load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; LCARTESIAN = False

begin
;=================================================;
; Open file
;=================================================;
; MESONH input file
  mnh_file="REUNI.1.00A20.004.nc"
  a = addfile(mnh_file, "r")

;==================================================;
; Open the workstation and choose colormap
; For paper quality plot do not use type ncgm or eps
; Use type= ps or pdf 
;==================================================;

  type = "x11"          ; could be png, ps , pdf ...
  wks = gsn_open_wks(type,"plt_Cloud")

; Choose the colormap you want 
; See http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
; for all the available color tables
  gsn_define_colormap(wks,"WhBlGrYeRe")   ; overwrite the .hluresfile color map

;=================================================;
; Get informations on variable sizes
;=================================================;
  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"UT") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2
; -2 : to remove non-physical values

;=================================================;
; Read the variables we will need
;=================================================;

; Read latitude and longitude fields
; -----------------------------------------
  lat2d = a->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lat2d@units="degrees_north"
  lon2d = a->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lon2d@units="degrees_east"

; Read water wapor mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RVT"))
      qv = a->RVT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
      qv = qv*1000.
      qv@units = "g/kg"   
      qv@description="Water vapor mixing ratio"
      qv@lat2d=lat2d
      qv@lon2d=lon2d
    end if

; Read cloud mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RCT"))
      qc = a->RCT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
      qc = qc*1000.
      qc@units = "g/kg"   
      qc@description="Cloud mixing ratio"
      qc@lat2d=lat2d
      qc@lon2d=lon2d
    end if

; Read rain mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RRT"))
      qr = a->RRT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
      qr = qr*1000.
      qr@units = "g/kg"   
      qr@description="Rain mixing ratio"
      qr@lat2d=lat2d
      qr@lon2d=lon2d
    end if

; Read ice mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RIT"))
      qi = a->RIT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
      qi = qi*1000.
      qi@units = "g/kg"   
      qi@description="Ice mixing ratio"
      qi@lat2d=lat2d
      qi@lon2d=lon2d
    end if

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

;=================================================;
; Set some other basic resources
;=================================================;

  resmap@mpDataBaseVersion      = "HighRes"  ; choose highres map data version (must be donwloaded)
  resmap@mpFillOn               = False

; Defining the corners for projection
;====================================
  resmap@mpLimitMode            = "Corners"
  resmap@mpLeftCornerLatF       = lat2d(0,0)
  resmap@mpLeftCornerLonF       = lon2d(0,0)
  resmap@mpRightCornerLatF     = lat2d(jmax-1,imax-1)
  resmap@mpRightCornerLonF     = lon2d(jmax-1,imax-1)


;=========================================================;
;  Loop over levels by step of 5
;=========================================================;
    do level = 0,19,5      

      display_level = level + 1

; Contour resources
;--------------------
      opts           =  resmap
      opts@cnFillOn  = True
      opts@cnLinesOn = False

      opts@tiMainString = "Model Level "+display_level


; Plot wapor mixing ratio if it exists
; -----------------------------------------
      if (isvar("qv"))
        plot=gsn_csm_contour_map(wks,qv(level,:,:),opts)
      end if

; Plot cloud mixing ratio if it exists
; -----------------------------------------
      if (isvar("qc"))
        plot=gsn_csm_contour_map(wks,qc(level,:,:),opts)
      end if

; Plot rain mixing ratio if it exists
; -----------------------------------------
      if (isvar("qr"))
        plot=gsn_csm_contour_map(wks,qr(level,:,:),opts)
      end if

; Plot ice mixing ratio if it exists
; -----------------------------------------
      if (isvar("qi"))
        plot=gsn_csm_contour_map(wks,qi(level,:,:),opts)
      end if


    end do      ; END OF LEVEL LOOP
;=================================================;

end
