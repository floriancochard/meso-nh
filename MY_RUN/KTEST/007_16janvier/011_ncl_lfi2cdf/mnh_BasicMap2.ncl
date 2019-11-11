;================================================;
; Always load these files in the same order
; -----------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "MESONHtools.ncl"   

begin
;=================================================;
; Open file 
;=================================================;
  mnh_file="16JAN.1.12B18.001.cdf"
  a = addfile(mnh_file, "r")

;=================================================;
; Get informations on variable sizes
;=================================================;
  mdims = getfilevardimsizes(a,"ZS") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2
  jmax=mdims(nd-2)-2

;=================================================;
; Read some variables frm mnh_file
;=================================================;
  ZS = a->ZS(1:jmax,1:imax)
  ZS@long_name="Orography"
  ZS@units="m"

  Z=a->ZHAT
  Y=a->YHAT(1:jmax)
  X=a->XHAT(1:imax)

;=================================================;
; Use MESONHtools procedure to set projection info
;=================================================;
  resmap                        = True          ; plot mods desired

  imax = a->IMAX
  jmax = a->JMAX
  lat2d = new((/jmax,imax/),"double")
  lat2d(:,:)=0.
  lon2d = new((/jmax,imax/),"double")
  lon2d(:,:)=0.
  icornerMap = new((/2,2/),"integer")
  icornerMap(:,:)=0
  MESONH_map_c (a, resmap, lat2d, lon2d, icornerMap)  

;=================================================;
; Map resources 
;=================================================;
  resmap@mpDataBaseVersion    	= "HighRes" ; choose highres map data version (must be donwloaded)
  resmap@mpGridAndLimbOn   	= True          ; turn on lat/lon lines
  resmap@mpGridLatSpacingF 	= 10.           ; spacing for lat lines
  resmap@mpGridLonSpacingF 	= 10.           ; spacing for lon lines
  resmap@mpGeophysicalLineColor = "Black"
  resmap@mpNationalLineColor    = "Black"
  resmap@mpGridLineColor        = "Black"
  resmap@mpLimbLineColor        = "Black"
  resmap@mpPerimLineColor       = "Black"

  resmap@mpGeophysicalLineThicknessF = 2.0

;=================================================;
; Plot Orography (ZS)  (filled contour)
;=================================================;

; Open workstation and define colormap 
; -------------------------------------
  wks  = gsn_open_wks("x11","plt_BasicMap2")   ; open a X11
  gsn_define_colormap(wks,"topo_15lev") ; Choose colormap

; General resources
; --------------------
  resmap@gsnMaximize               = True      ;  to maximize size
  resmap@gsnSpreadColors    	= True      ; use full range of colormap

; Contour resources
; --------------------
  resmap@cnFillOn  		=  True        ; turn on color fill
  resmap@cnLinesOn 		= False          ; turn off contour lines
  resmap@cnLevelSelectionMode = "ExplicitLevels"
  cnLevels  = (/1,25,50,75,100,125,150,175,200,225,250,275,300/)
  resmap@cnLevels    = cnLevels

  resmap@lbOrientation  = "Vertical"

; Plot 
; -----
  plot = gsn_csm_contour_map(wks,ZS,resmap)

end


