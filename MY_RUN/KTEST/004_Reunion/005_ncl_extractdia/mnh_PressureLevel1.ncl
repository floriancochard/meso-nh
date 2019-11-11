load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "MESONHtools.ncl"

;======================================================================;
; Example script to interpolate from model to pressure level and 
; to plot some fields from MESONH netcdf file generated with extractdia
;======================================================================;

begin

;=======================================================;
; Pressure levels that we want the data interpolated to
;=======================================================;
  pressure_levels = (/ 850., 700., 500., 300./)
  nlevels         = dimsizes(pressure_levels)  ; number of pressure levels

;=================================================;
; Open file 
; MESONH netcdf input file generated with extractdia
;=================================================;
  mnh_file="REUNI.1.00A20.004dgKCL.nc"
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
; Read the variables we need
;=================================================;
    tk = a->THT(0,:,:,:)	          ; potential temperature
    tk@description="Potential Temperature"

    u = a->UT(0,:,:,:)               ; u

    v = a->VT(0,:,:,:)               ; v

    p = a->PABST(0,:,:,:)            ; pressure
    p@description="Pressure"

    rv = a->RVT(0,:,:,:)*1000        ; water vapor mixing ratio
    rv@description="Water vapor mixing ratio"

    ZS = a->ZSBIS(0,:,:) 	     ; terrain

    z = a->VLEV(:,:,:)		     ; heigth 
    z@description="Height"


;=================================================;
; Horizontal interpolation on pressure levels
;=================================================;
    do level = 0,nlevels-1                 ; LOOP OVER LEVELS

      pressure = pressure_levels(level)

      tk_plane = wrf_user_intrp3d(tk,p,"h",pressure,0.,False)
      z_plane  = wrf_user_intrp3d( z,p,"h",pressure,0.,False)
      u_plane  = wrf_user_intrp3d( u,p,"h",pressure,0.,False)
      v_plane  = wrf_user_intrp3d( v,p,"h",pressure,0.,False)
      rv_plane = wrf_user_intrp3d( rv,p,"h",pressure,0.,False)

      spd_plane     = (u_plane*u_plane + v_plane*v_plane)^(0.5) ; m/sec
      spd_plane@description = "Wind Speed"
      spd_plane@units = "m/s"

      u_plane@units = "m/s"
      v_plane@units = "m/s"


;==================================
; Create plots using wrf functions
;==================================
      res=True
      res@MainTitle                   = "MESONH"
      res@Footer = False


; Plotting options for Tk                
;------------------------
        opts = res                          
        opts@cnLineColor = "Red"
;        opts@ContourParameters = (/ 0.1 /)
        opts@cnInfoLabelOrthogonalPosF = 0.07  ; offset second label information
        opts@gsnContourLineThicknessesScale = 2.0
        contour_tk = wrf_contour(a,wks,tk_plane,opts)
        delete(opts)


; Plotting options for  RV          
;--------------------------
        opts                          = res                          

        opts@cnFillOn                 = True  
        opts@cnMissingValFillColor    = "Red"
        opts@cnMissingValFillPattern  = 8

        opts@pmLabelBarOrthogonalPosF = -0.1

        contour_rv = wrf_contour(a,wks,rv_plane,opts)
        delete(opts)


; Plotting options for Wind Speed                
;--------------------------------
        opts                                = res                          

        opts@gsnContourLineThicknessesScale = 3.0

        opts@cnLineColor                    = "MediumSeaGreen"
        opts@cnInfoLabelOrthogonalPosF      = 0.07  ; offset second label information

        contour_spd = wrf_contour(a,wks,spd_plane,opts)
        delete(opts)


; Plotting options for Wind Vectors                 
;----------------------------------
        opts            = res          

        opts@FieldTitle = "Wind"   ; overwrite Field Title
        opts@NumVectors = 47       ; wind barb density

        vector = wrf_vector(a,wks,u_plane,v_plane,opts)
        delete(opts)


; Plotting options for Geopotential Heigh
;----------------------------------------
        opts_z                                = res                          
        opts_z@gsnContourLineThicknessesScale = 2.0
        opts_z@cnLineColor                    = "Blue"


;============
; MAKE PLOTS                                       
;============
        pltres = True

; Resources for map
;---------------------
        resmap = True
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
