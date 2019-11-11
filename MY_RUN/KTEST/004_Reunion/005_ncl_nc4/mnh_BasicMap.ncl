load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;=============================================================;
;  Example of script to plot ZS field from MESONH netcdf4 file
;  Could be use with any files generated by MESONH if LCARTESIAN=False
;=============================================================;
begin
;=============================================================;
; Open file 
;=============================================================;

  mnh_file="REUNI.1.00A20.004.nc"
  a = addfile(mnh_file, "r")

;=================================================;
; Get information on variable sizes
;=================================================;

  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"ZS") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext  ; -2*jphext : to remove non-physical values
  jmax=mdims(nd-2)-2*jphext  ; -2*jphext : to remove non-physical values

;=================================================;
; Read the variables 
;=================================================;

; Note: do not read first and last value which are non physical values
; -------------------------------------

  lat2d=a->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  lon2d=a->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)

  ZS = a->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1)
  ZS@long_name="Orography"
  ZS@units="m"

; Associate coordinates to variable
; Necessary to plot the field at the right geographical position 
; -------------------------------------

  ZS@lat2d = lat2d
  ZS@lon2d = lon2d

; Read projection parameters 
; -------------------------------------

  RPK = a->RPK
  BETA =a->BETA
  LON0= a->LON0
  
;=================================================;
; CREATE PLOT : Orography (ZS, filled contour)
;=================================================;
; Open workstation and define colormap
; -------------------------------------

  type  = "x11" ; open a x11 window
                ; change type to png, ps or pdf to get the plot into a file
  wks  = gsn_open_wks(type,"plt_BasicMap")   ; 
  gsn_define_colormap(wks,"topo_15lev") ; Choose colormap


;==========================================
; Map ressources
;==========================================
  resmap                        = True         

; If there is an error on HighRes, it means that you don't have the HighRes data
; You need to download them or the change HighRes by MediumRes
; See https://www.ncl.ucar.edu/Document/Graphics/rangs.shtml for info
; -------------------------------------
  resmap@mpDataBaseVersion    	= "HighRes"	; choose highres map data version (must be donwloaded)
  resmap@mpGridAndLimbOn   	= True             ; turn on lat/lon lines
  resmap@mpGridLatSpacingF 	= 10.              ; spacing for lat lines
  resmap@mpGridLonSpacingF 	= 10.              ; spacing for lon lines

  resmap@pmTickMarkDisplayMode = "Always"         ; turn on tickmarks

;=================================================;
; Set map projection ressources using projection parameters
;=================================================;
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

;==========================================
; Defining the corners for projection otherwise the plot will be global
;==========================================
  resmap@mpLimitMode           = "Corners"
  resmap@mpLeftCornerLatF      = lat2d(0,0)
  resmap@mpLeftCornerLonF      = lon2d(0,0)
  resmap@mpRightCornerLatF     = lat2d(jmax-1,imax-1)
  resmap@mpRightCornerLonF     = lon2d(jmax-1,imax-1)


 print ("Corner (0,0); Lat="+resmap@mpLeftCornerLatF+ \
                    ", Lon="+resmap@mpLeftCornerLonF)
 print ("Oppos corner; Lat="+resmap@mpRightCornerLatF+ \
                     ", Lon= "+resmap@mpRightCornerLonF)


;==========================================
; Create ZS plot  (contour)
;==========================================
; ---------------------------
; General ressources
; ---------------------------
  resmap@gsnMaximize            = True          ; maximize size
  resmap@gsnSpreadColors    	= True       	; use full range of colormap

; ---------------------------
; Contour ressources
; ---------------------------
  resmap@cnFillOn  		= True          ; turn on color fill
  resmap@cnLinesOn 		= False         ; turn off contour lines
  resmap@cnLevelSelectionMode 	= "ManualLevels"; Manually set contour levels.
  resmap@cnMinLevelValF       	=  5.
;  resmap@cnMaxLevelValF       	= 3200.
;  resmap@cnLevelSpacingF      	=  200.

;=================================================;
; PLOT ZS 
;=================================================;
  plot_zs = gsn_csm_contour_map(wks,ZS,resmap)

end

