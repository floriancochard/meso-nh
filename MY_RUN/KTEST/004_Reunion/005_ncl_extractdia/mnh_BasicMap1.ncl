;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "MESONHtools.ncl"   

; ================================================;
begin
;=================================================;
; Open MESONH file 
;=================================================;
  mnh_file="REUNI.1.00A20.004dgKCL.nc"
  a = addfile(mnh_file, "r")

;=================================================;
; Get informations on variable sizes
;=================================================;
  mdims = getfilevardimsizes(a,"ZSBIS") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  IMAX=mdims(nd-1)
  JMAX=mdims(nd-2)

;=================================================;
; Read some variables
;=================================================;
  ZS = a->ZSBIS(0,:,:)
  ZS@long_name="Orography"
  ZS@units="m"

;==============================================
; Use MESONHtools procedure to set projection 
;==============================================
  resmap     = True

  lat2d = new((/JMAX,IMAX/),"double")
  lat2d(:,:)=0.
  lon2d = new((/JMAX,IMAX/),"double")
  lon2d(:,:)=0.
  icornerMap = new((/2,2/),"integer")
  icornerMap(:,:)=0
  MESONH_map_c (a, resmap, lat2d, lon2d, icornerMap)  

;=================================================;
; More map resources
;=================================================;
  resmap@mpDataBaseVersion      = "HighRes" ; highres version (must be donwloaded)
  resmap@mpGridAndLimbOn        = True          ; turn on lat/lon lines
  resmap@mpGridLatSpacingF      = 0.5           ; spacing for lat lines
  resmap@mpGridLonSpacingF      = 0.5           ; spacing for lon lines
  resmap@mpGeophysicalLineColor = "Black"
  resmap@mpNationalLineColor    = "Black"
  resmap@mpGridLineColor        = "Black"
  resmap@mpLimbLineColor        = "Black"
  resmap@mpPerimLineColor       = "Black"

;=================================================;
; PLOT Orography (filled contour)
;=================================================;

; Open workstation and define colormap
; -------------------------------------
  wks  = gsn_open_wks("X11","plt_BasicMap1")   ; open a ncgm file
  gsn_define_colormap(wks,"topo_15lev") ; Choose colormap

; General resources
; --------------------
  resmap@gsnMaximize               = True      ;  to maximize size
  resmap@gsnSpreadColors           = True      ; use full range of colormap

; Contour resources
;-------------------
  resmap@cnFillOn  		=  True        ; turn on color fill
  resmap@cnLinesOn 		= False          ; turn off contour lines
  resmap@cnLevelSelectionMode = "ExplicitLevels"
  cnLevels  = (/1,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000/)
  resmap@cnLevels    = cnLevels

  resmap@lbLabelAngleF    =    45                 ; angle labels



;=================================================;
; PLOT ZS 
;=================================================;
  plot = gsn_csm_contour_map(wks,ZS,resmap)

end


