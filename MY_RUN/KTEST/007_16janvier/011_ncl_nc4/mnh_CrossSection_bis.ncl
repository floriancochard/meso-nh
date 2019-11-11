load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "MESONHtools.ncl"

begin

;=================================================;
; Open file 
;=================================================;
  mnh_file="16JAN.1.12B18.001.nc"
  a = addfile(mnh_file, "r")

;=================================================;
; We generate plots, but what kind do we prefer?
; For paper quality plot do not use ncgm or eps
;=================================================;
  type = "x11"

;=================================================;
; Open the workstation and choose colormap
;=================================================;
  wks = gsn_open_wks(type,"mnh_CrossSection_bis")
;  gsn_define_colormap(wks,"rainbow")

;=================================================;
; Get informations on variable sizes
;=================================================;
  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"UT") ; get some dimension sizes for the file
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2

;=================================================;
; Read  some variables from mnh_file
;=================================================;
    lat2d = a->LAT(jphext:jmax+jphext-1,jphext:imax+jphext-1)
    lat2d@units="degrees_north"
    lon2d = a->LON(jphext:jmax+jphext-1,jphext:imax+jphext-1)
    lon2d@units="degrees_east"

    ZS =a->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1)              
    ZS@lat2d = lat2d
    ZS@lon2d = lon2d

    rv = a->RVT(1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
    rv = rv*1000.
    rv@description="RVT"
    rv@lat2d = lat2d
    rv@lon2d = lon2d

;    zh= a->ZHAT(1:kmax+1)
    zh= a->ZHAT(0:kmax)

; Read parameter for projection

   RPK  = a->RPK
   BETA = a->BETA
   LON0 = a->LON0

;=================================================;
; Take care of the Arakawa C grid
;=================================================;
; Unstagger zh (from grid 4 to 1)
;----------------------------------;
    nzh=new(kmax,double)
    do k=0,kmax-1
     nzh(k)=(zh(k)+zh(k+1))/2.
    end do

;=================================================;
; Needed for the wrf_intrp3d function
;=================================================;
; Create Z3D 
;----------------------------------;
    z=new(dimsizes(rv),double)
    zcoef=1.-ZS/nzh(kmax-1)

    do i=0,imax-1
    do j=0,jmax-1
       z(:,j,i) = nzh*zcoef(j,i)+ZS(j,i)
    end do
    end do

; get height info for labels
;----------------------------------;
      zmin = 0.
      zmax = 6.         ; We are only interested in the first 6km
      nz   = floattoint(zmax + 1)

;=================================================;
; Set map projection
;=================================================;

; Resources for the map projection
;-----------------------------------------------;
  resmap=True
;-----------------------------------------------;
; Get global attributes to set map projection
;-----------------------------------------------;

   if (RPK.gt.0)
  ;   Lambert  projection from north pole
; ---------------------------
   resmap@mpProjection          = "LambertConformal"     ; projection
   pole = 1 ; with pole=1 for north hemisphere and -1 for south hemisphere
   resmap@mpLambertParallel1F   =  pole*asin(RPK)*180/3.14 ; 
   resmap@mpLambertParallel2F   = resmap@mpLambertParallel1F  ;  Parallel1 = PArallel2
   resmap@mpLambertMeridianF    =  LON0      ; get value from globla attribute
  end if

  if (RPK.lt.0)
  ;   Lambert projection from south pole
; ---------------------------
   resmap@mpProjection          = "LambertConformal"     ; projection
   pole = -1 ; with pole=1 for north hemisphere and -1 for south hemisphere
   resmap@mpLambertParallel1F   =  pole*asin(RPK)*180/3.14 ; 
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

;==========================================
; Needed for regional native projection
;==========================================
  resmap@gsnAddCyclic         = False              ; regional data

;=================================================;
; Define Other Resources for map ploting
;=================================================;

  resmap@gsnFrame               = False
  resmap@gsnMaximize            = True

; Map resources
;-----------------------------------------------;
  resmap@mpFillOn               = False
  resmap@mpGridAndLimbOn        = True            ; turn on lat/lon lines
  resmap@mpGridLatSpacingF      = 10.             ; spacing for lat lines
  resmap@mpGridLonSpacingF      = 10.             ; spacing for lon lines
  resmap@mpDataBaseVersion      = "HighRes"     ; choose highres map data version (must be donwloaded)
  resmap@pmTickMarkDisplayMode = "Always"         ; turn on tickmarks
  resmap@mpOutlineBoundarySets = "AllBoundaries" ; state boundaries
  resmap@mpPerimDrawOrder      = "PostDraw"       ; force map perim


; Contour definitions
;-----------------------------------------------;
  resmap@cnFillOn               = True
  resmap@cnLinesOn              = False
;  resmap@cnLevelSelectionMode   = "ExplicitLevels"
;  cnLevels                      = (/1,25,50,75,100,125,150,175,200,225,250,275,300/)
;  resmap@cnLevels               = cnLevels

  resmap@lbOrientation          = "Vertical"

  
  FirstPlot = True

;=================================================;
; Planes definition and interpolation
;=================================================;
    do ip = 1, 3      	; we are doing 3 plots
			; all with the pivot point (plane) in the center of the domain
			; at angles 0, 45 and 90
 ;  
 ;                   |
 ;       angle=0 is  |
 ;                   |
 ; 

; Define the plane where the cross-section will be done
; ------------------------------------------------------
        plane = new(2,float)
        plane = (/ (imax/2)+10, (jmax)/2 /)    
; Resource for the interpolation
; --------------------------------
        opts = False

        if(ip .eq. 1) then
          angle = 90.
          X_plane = wrf_user_intrp2d(lon2d,plane,angle,opts)
          X_desc = "longitude"
        end if
        if(ip .eq. 2) then
          angle = 0.
          X_plane = wrf_user_intrp2d(lat2d,plane,angle,opts)
          X_desc = "latitude"
        end if
        if(ip .eq. 3) then
          angle = -35.
          X_plane = wrf_user_intrp2d(lon2d,plane,angle,opts)
          X_desc = "longitude"
          
        end if

; Interpolate fields onto the defined plane
; --------------------------------------------
        rv_plane = wrf_user_intrp3d(rv,z,"v",plane,angle,opts)

; Find the index where 6km is - only need to do this once
; --------------------------------------------
        if ( FirstPlot ) then
          zz = wrf_user_intrp3d(z,z,"v",plane,angle,opts)
          b = ind(zz(:,0) .gt. zmax*1000. )
          zmax_pos = b(0) - 1
          if ( abs(zz(zmax_pos,0)-zmax*1000.) .lt. abs(zz(zmax_pos+1,0)-zmax*1000.) ) then
            zspan = b(0) - 1
          else
            zspan = b(0)
          end if
          delete(zz)
          delete(b)
          FirstPlot = False
        end if

; Define X-axis labels
;----------------------
      dimsX = dimsizes(X_plane)
      xmin  = X_plane(0)
      xmax  = X_plane(dimsX(0)-1)
      xspan = dimsX(0)-1
      nx    = 3

; Options for XY Plots
;-----------------------------------------
        opts_xy                         = True
        opts_xy@tiXAxisString           = X_desc
        opts_xy@tiYAxisString           = "Height (km)"
        opts_xy@tiXAxisFontHeightF      = 0.020
        opts_xy@tiYAxisFontHeightF      = 0.020

; Resources to plot the topography (missing values in fields)
;-----------------------------------------
        opts_xy@cnMissingValPerimOn     = True
        opts_xy@cnMissingValFillColor   = "red"
        opts_xy@cnMissingValFillPattern = 8
    
; Tickmarks/labels resources
;------------------------------
        opts_xy@lbOrientation           = "Vertical"

        opts_xy@tmXTOn                  = False
        opts_xy@tmYROn                  = False
        opts_xy@tmXBMode                = "Explicit"
        opts_xy@tmXBValues              = fspan(0,xspan,nx)              ; Create tick marks
        opts_xy@tmXBLabels              = sprintf("%.2f",fspan(xmin,xmax,nx))  ; Create labels
        opts_xy@tmXBLabelFontHeightF    = 0.015
        opts_xy@tmYLMode                = "Explicit"
        opts_xy@tmYLValues              = fspan(0,zspan,nz)                    ; Create tick marks
        opts_xy@tmYLLabels              = sprintf("%.1f",fspan(zmin,zmax,nz))  ; Create labels
        opts_xy@tmXBMajorLengthF        = 0.02
        opts_xy@tmYLMajorLengthF        = 0.02
        opts_xy@tmYLLabelFontHeightF    = 0.015


; Plotting options for RV
;--------------------------
        opts_rv                         = opts_xy
        opts_rv@cnFillOn                = True
        opts_rv@gsnFrame = False
        opts_rv@gsnDraw  = False
        opts_rv@cnLevelSelectionMode = "ManualLevels"
        opts_rv@cnMinLevelValF        =  0.
        opts_rv@cnMaxLevelValF        = 10.
        opts_rv@cnLevelSpacingF       =  0.5

;======================
; MAKE PLOTS         
;======================

; Horizontal plot with topography
;-----------------------------------------------------------------------
        plot=gsn_csm_contour_map(wks,ZS,resmap)

; Draw line that shows the position of the cross section
;-----------------------------------------------------------------------
        lat_plane = wrf_user_intrp2d(lat2d,plane,angle,opts)
        lon_plane = wrf_user_intrp2d(lon2d,plane,angle,opts)

        lnres = True
        lnres@gsLineThicknessF = 3.0
        lnres@gsLineColor = "Red"
        do ii = 0,dimsX(0)-2
            gsn_polyline(wks,plot,(/lon_plane(ii),lon_plane(ii+1)/),(/lat_plane(ii),lat_plane(ii+1)/),lnres)
        end do
        frame(wks)
        delete(lon_plane)
        delete(lat_plane)

; Vertical cross-section
;--------------------------



        contour_rv = gsn_csm_contour(wks,rv_plane(0:zmax_pos,:),opts_rv)

        draw(contour_rv)
        frame(wks)


; Delete options and fields, so we don't have carry over
;-------------------------------------------------------
        delete(opts_xy)
        delete(opts_rv)
        delete(rv_plane)
        delete(X_plane)

end do  ; make next cross section

end
