
;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "MESONHtools.ncl"

begin
;=================================================;
; Open file
;=================================================;
; MESONH input file
  mnh_file="16JAN.1.12B18.001.cdf"
  a = addfile(mnh_file, "r")

;==================================================;
; Open the workstation and choose colormap
; For paper quality plot do not use type ncgm or eps
; Use type= ps or pdf or x11 (for debugging)
;==================================================;

  type = "x11"
  wks = gsn_open_wks(type,"plt_Cloud2")
  gsn_define_colormap(wks,"WhBlGrYeRe")   ; overwrite the .hluresfile color map

;=================================================;
; Get informations on variable sizes
;=================================================;

  mdims = getfilevardimsizes(a,"RVT") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2
  jmax=mdims(nd-2)-2
  kmax=mdims(nd-3)-2
; -2 : to remove non-physical values

;=================================================;
; Read the variables we will need
;=================================================;

; Read water wapor mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RVT"))
      qv = a->RVT(1:kmax,1:jmax,1:imax)
      qv = qv*1000.
      qv@units = "g/kg"   
      qv@description="Water vapor mixing ratio"
    end if

; Read cloud mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RCT"))
      qc = a->RCT(1:kmax,1:jmax,1:imax)
      qc = qc*1000.
      qc@units = "g/kg"   
      qc@description="Cloud mixing ratio"
    end if

; Read rain mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RRT"))
      qr = a->RRT(1:kmax,1:jmax,1:imax)
      qr = qr*1000.
      qr@units = "g/kg"   
      qr@description="Rain mixing ratio"
    end if

; Read ice mixing ratio if it exists
; -----------------------------------------
    if(isfilevar(a,"RIT"))
      qi = a->RIT(1:kmax,1:jmax,1:imax)
      qi = qi*1000.
      qi@units = "g/kg"   
      qi@description="Ice mixing ratio"
    end if

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
;-----------------------------------------------;
  MESONH_map_c (a, resmap, lat2d, lon2d, icornerMap)

;=================================================;
; Set some other basic resources
;=================================================;

  resmap@mpDataBaseVersion      = "HighRes"  ; choose highres map data version (must be donwloaded)
  resmap@mpFillOn               = False

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
