load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "MESONHtools.ncl"

;=================================================;
; Example of script to interpolate from model to pressure level 
; and to plot some fields from MESONH netcdf file generated with lfi2cdf
;=================================================;

begin

;=======================================================;
; Pressure levels that we want the data interpolated to
;=======================================================;
  pressure_levels = (/ 850., 700., 500., 300./)  
  nlevels         = dimsizes(pressure_levels)  ; number of pressure levels

;=================================================;
; Open file 
; MESONH netcdf input file generated with lfi2cdf  
;=================================================;

  mnh_file="16JAN.1.12B18.001.cdf"
  a = addfile(mnh_file, "r")

;=================================================;
; Open the workstation and choose colormap
; For paper quality plot do not use type ncgm or eps
; Use ps or pdf and x11 (for debugging)
;==================================================;

  type = "x11"
  wks = gsn_open_wks(type,"plt_PressureLevel2")
  gsn_define_colormap(wks,"spread_15lev")

;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  mdims = getfilevardimsizes(a,"UT") ; get dimension sizes
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2
  jmax=mdims(nd-2)-2
  kmax=mdims(nd-3)-2

;=================================================;
; Read the variables we need
;=================================================;

  tk = a->THT(1:kmax,1:jmax,1:imax)	     ; potential temperature
  tk@description="Potential Temperature"

  u = a->UT(1:kmax,1:jmax,1:imax+1)         ; u
  v = a->VT(1:kmax,1:jmax+1,1:imax)         ; v

; Unstagger u,v
  u!0="Z"
  u!1="Y"
  u!2="X"
  u_unst = wrf_user_unstagger(u,"X")

  v!0="Z"
  v!1="Y"
  v!2="X"
  v_unst = wrf_user_unstagger(v,"Y")

  p = a->PABST(1:kmax,1:jmax,1:imax)        ; pressure
  p@description="Pressure"

  rv = a->RVT(1:kmax,1:jmax,1:imax)*1000    ; water vapor mixing ratio
  rv@description="Water vapor mixing ratio"

  zs = a->ZS(1:jmax,1:imax) 	            ; terrain
  zh = a->ZHAT(1:kmax+1)	            ; heigth without terrain

;==========================================
; Create 3D Z to be used for interpolation
;==========================================

  z=new(dimsizes(tk),double)
  z@description="Height"

; Unstag zh
;----------
   nzh=new(kmax,double)
   do k=0,kmax-1
   nzh(k)=(zh(k)+zh(k+1))/2.
   end do
   zcoef=1.-zs/nzh(kmax-1)

    do i=0,imax-1
    do j=0,jmax-1
       z(:,j,i) = nzh*zcoef(j,i)+zs(j,i)
    end do
    end do

;=================================================;
; Horizontal interpolation on pressure levels
;=================================================;


    do level = 0,nlevels-1                 ; LOOP OVER LEVELS

      pressure = pressure_levels(level)

      tk_plane = wrf_user_intrp3d(tk,p,"h",pressure,0.,False)
      z_plane  = wrf_user_intrp3d( z,p,"h",pressure,0.,False)
      u_plane  = wrf_user_intrp3d( u_unst,p,"h",pressure,0.,False)
      v_plane  = wrf_user_intrp3d( v_unst,p,"h",pressure,0.,False)
      rv_plane  = wrf_user_intrp3d( rv,p,"h",pressure,0.,False)

      spd_plane = (u_plane*u_plane + v_plane*v_plane)^(0.5) ; m/sec
      spd_plane@description = "Wind Speed"
      spd_plane@units = "m/s"

      u_plane@units = "m/s"
      v_plane@units = "m/s"

;==================================
; Create plots usinf wrf functions
;==================================
      res=True
      res@MainTitle                   = "MESONH"
      res@Footer = False

; Plotting options for Tk                
;--------------------------
      opts =  res                         
      opts@cnLineColor = "Red"
      opts@cnInfoLabelOrthogonalPosF = 0.07 ; offset second label information
      opts@gsnContourLineThicknessesScale = 2.0

      contour_tk = wrf_contour(a,wks,tk_plane,opts)
      delete(opts)

; Plotting options for  RV          
;--------------------------
  opts                          = res                          
  opts@cnFillOn                 = True  
  opts@pmLabelBarOrthogonalPosF = -0.1
  opts@cnMissingValFillColor    = "Red"
  opts@cnMissingValFillPattern  = 8

  contour_rv = wrf_contour(a,wks,rv_plane,opts)
  delete(opts)

; Plotting options for Wind Speed                
;--------------------------
  opts = res                          
  opts@cnLineColor = "MediumSeaGreen"
  opts@cnInfoLabelOrthogonalPosF = 0.07  ; offset second label information
  opts@gsnContourLineThicknessesScale = 3.0

  contour_spd = wrf_contour(a,wks,spd_plane,opts)
  delete(opts)


; Plotting options for Wind Vectors                 
;-----------------------------------
  opts            = res          
  opts@FieldTitle = "Wind"   ; overwrite Field Title
  opts@NumVectors = 47       ; wind barb density

  vector = wrf_vector(a,wks,u_plane,v_plane,opts)
  delete(opts)

; Plotting options for Geopotential Heigh
;-----------------------------------------
  opts_z                                = res                          
  opts_z@cnLineColor                    = "Blue"
  opts_z@gsnContourLineThicknessesScale = 2.0

;=============
; MAKE PLOTS                                       
;=============

  pltres = True

; Resources for map
;-------------------
  resmap                             = True
  resmap@mpGeophysicalLineColor      = "Black"
  resmap@mpNationalLineColor         = "Black"
  resmap@mpGridLineColor             = "Black"
  resmap@mpLimbLineColor             = "Black"
  resmap@mpPerimLineColor            = "Black"
  resmap@mpGeophysicalLineThicknessF = 2.0
  resmap@mpGridLineThicknessF        = 2.0
  resmap@mpLimbLineThicknessF        = 2.0
  resmap@mpNationalLineThicknessF    = 2.0
  resmap@mpGridAndLimbOn             = False            ; turn off lat/lon lines
  resmap@mpDataBaseVersion           = "HighRes"     ; choose highres map data


; Overlay rv , tk , height and wind barbs
;------------------------------------------
  if ( pressure .eq. 850 ) then   ; plot temp, rv, height, wind barbs

   contour_height = wrf_contour(a,wks,z_plane,opts_z)
   plot = mnh_map_overlays(a,wks,(/contour_rv,contour_tk,contour_height, \
                                    vector/),pltres,resmap)
  end if

; Overlay  tk , height and wind barbs
;------------------------------------------
  if ( pressure .eq. 700 ) then   ; plot temp, height, wind barbs

   contour_height = wrf_contour(a,wks, z_plane,opts_z)
   plot = mnh_map_overlays(a,wks,(/contour_tk,contour_height, \
                                    vector/),pltres,resmap)
  end if

; Overlay  tk , height and wind barbs
;------------------------------------------
  if ( pressure .eq. 500 ) then   ; plot temp, height, wind barbs

   contour_height = wrf_contour(a,wks, z_plane,opts_z)
   plot = mnh_map_overlays(a,wks,(/contour_tk,contour_height, \
                                    vector/),pltres,resmap)
  end if

; Overlay  wind speed , height and wind barbs
;------------------------------------------
  if ( pressure .eq. 300 ) then   ; plot windspeed, height, wind barbs

   contour_height = wrf_contour(a,wks, z_plane,opts_z)
   plot = mnh_map_overlays(a,wks,(/contour_spd,contour_height, \
                                    vector/),pltres,resmap)
  end if

  delete(opts_z)

end do      ; END OF LEVEL LOOP


end
