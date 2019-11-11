load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;-------------------------------------------------------------
;contains:
; procedure MESONH_map_c
;function mnh_map(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical)   
;function mnh_map_overlays(in_file[1]:file,wks:graphic,plots[*]:graphic, \
;                                  opt_arg[1]:logical,opt_mp[1]:logical)
;function MESONH_pinter( pfield:numeric, loc_param:numeric, ppabs:numeric )
;-------------------------------------------------------------

;==============================================================
; J.-P. CHABOUREAU 
; This is a driver that selects the appropriate 
; mapping function based upon the file variables RPK, BETA, LATOR, LONOR
;
;
; Sample usage:
;             a = addfile("...", r")
;             IMAX = a->IMAX
;             JMAX = a->JMAX
;             lat2d = new((/JMAX,IMAX/),"double")
;             lat2d(:,:)=0.
;             lon2d = new((/JMAX,IMAX/),"double")
;             lon2d(:,:)=0.
;             icorners = new((/2,2/),"integer")
;             icorners(:,:)=0
;             res  = True
;             MESONH_map_c (a, res, lat2d, lon2d, icorners)
;
;
undef("MESONH_map_c")
;==============================================================
procedure MESONH_map_c (in_file:file, res:logical, plat, plon, icorner)
;==============================================================
;local rank, dimll, nlat, mlon, lat, lon
local rank, dimll, nlat, mlon
begin

; Check if the variable RPK is in the file
; ----------------------------------------
if(isfilevar(in_file,"RPK")) then

; Read projection parameters
; -------------------------
   ZRPK   = in_file->RPK
   ZLATOR = in_file->LATOR
   ZLONOR = in_file->LONOR
   ZBETA = in_file->BETA
   ZLAT0 = in_file->LAT0
   ZLON0 = in_file->LON0

; Case netcdf from lfi2cdf
; -------------------------

   if(isfilevar(in_file,"IMAX"))
       XHAT=in_file->XHAT
       YHAT=in_file->YHAT
       IMAX= dimsizes(XHAT)-2
       JMAX= dimsizes(YHAT)-2
       zdx=XHAT(2)-XHAT(1)
       zdy=YHAT(2)-YHAT(1)

; unstagger
       do ji=0,IMAX-1
         XHAT(ji)=XHAT(ji)+zdx*1.5
       end do
       do jj=0,JMAX-1
         YHAT(jj)=YHAT(jj)+zdy*1.5
       end do

   else

; Case netcdf from extractdia
; ---------------------------
       XHAT=in_file->W_E_direction
       YHAT=in_file->S_N_direction
       IMAX= dimsizes(XHAT)
       JMAX= dimsizes(YHAT)
       zdx=XHAT(2)-XHAT(1)
       zdy=YHAT(2)-YHAT(1)

   end if

  print ("LATOR="+ZLATOR+" - LONOR="+ZLONOR)
  print ("ZLAT0="+ZLAT0+" - ZLON0="+ZLON0)
  print ("ZDX="+zdx+" - RPK="+ZRPK+" - BETA="+ZBETA)
  print ("IMAX="+IMAX+" - JMAX="+JMAX)

  if (ZRPK.gt.0) 
  ;   Stereographic projection
; ---------------------------
    res@mpProjection = "Stereographic"
    res@mpCenterLonF           = ZLON0
    res@mpCenterRotF = ZBETA
    res@mpCenterLatF           = 90. 
  end if

  if (ZRPK.lt.0) 
  ;   Stereographic projection
; ---------------------------
    res@mpProjection = "Stereographic"
    res@mpCenterLonF           = ZLON0
    res@mpCenterRotF = ZBETA
    res@mpCenterLatF           = -90. 
  end if

  if (ZRPK.eq.0) then
  ;   Mercator projection
; ---------------------------
    res@mpProjection = "Mercator"
  end if

 print("Map projection="+res@mpProjection)

else 
  print ("MESONH_map_c: Error no RPK variable in input file")
end if

;=================================================;
; calculate 2D lat and lon 
; based on src/mesonh_MOD/mode_gridproj.f90
;=================================================;

; Constants  
; -----------
  if(isfilevar(in_file,"IMAX"))
  XRADIUS=6371229.0d ; Earth radius (meters)
  else
  XRADIUS=6371.2290d ; Earth radius (km)
  end if
  XPI=2.0d*asin(1.)    ; Pi
  ZRDSDG= XPI/180.0d          ; Radian to Degree conversion factor
  ZXBM0 = 0.0d
  ZYBM0 = 0.0d

;=================================================;
  if (ZRPK.eq.0) then
; MERCATOR
;=================================================;
    XBETA=0.
    XLAT0=0.   ; map reference latitude  (degrees)
    ZXBM0 = 0.
    ZYBM0 = 0.
    ZCGAM    = cos(-ZRDSDG*XBETA)
    ZSGAM    = sin(-ZRDSDG*XBETA)
    ZRACLAT0 = XRADIUS*cos(ZRDSDG*ZLAT0)
    do ji=0,IMAX-1
      jj=0
      ZXMI0 = XHAT(ji)-ZXBM0
      ZYMI0 = YHAT(jj)-ZYBM0
      zlon  = (ZXMI0*ZCGAM+ZYMI0*ZSGAM)/(ZRACLAT0*ZRDSDG)+ZLONOR
      do jj=0,JMAX-1
        plon(jj,ji)=zlon
      end do
    end do
    do jj=0,JMAX-1
      ji=0
      ZXMI0 = XHAT(ji)-ZXBM0
      ZYMI0 = YHAT(jj)-ZYBM0
      ZT1  = log(tan(XPI/4.+ZLATOR*ZRDSDG/2.))
      ZT2  = (-ZXMI0*ZSGAM+ZYMI0*ZCGAM)/ZRACLAT0
      zlat  = (-XPI/2.+2.*atan(exp(ZT1+ZT2)))/ZRDSDG
      do ji=0,IMAX-1
        plat(jj,ji)=zlat
      end do
    end do
    
;=================================================;
  else
; STEREOGRAPHIC PROJECTION
;=================================================;
    ZCLAT0  = cos(ZRDSDG*ZLAT0)
    ZSLAT0  = sin(ZRDSDG*ZLAT0)
    ZCLATOR = cos(ZRDSDG*ZLATOR)
    ZSLATOR = sin(ZRDSDG*ZLATOR)
    ZRO0    = (XRADIUS/ZRPK)*(abs(ZCLAT0))^(1.-ZRPK) * \
              ((1.+ZSLAT0)*abs(ZCLATOR)/(1.+ZSLATOR))^ZRPK
    ZGA0    = (ZRPK*(ZLONOR-ZLON0)-ZBETA)*ZRDSDG
    ZXP     = ZXBM0-ZRO0*sin(ZGA0)
    ZYP     = ZYBM0+ZRO0*cos(ZGA0)
    do ji=0,IMAX-1
      do jj=0,JMAX-1
        ZATA = atan2( -(ZXP-XHAT(ji)) , (ZYP-YHAT(jj)) )/ZRDSDG
        zlon  = (ZBETA+ZATA)/ZRPK+ZLON0
        plon(jj,ji)=zlon
        ZRO2 = (XHAT(ji)-ZXP)^2+(YHAT(jj)-ZYP)^2
        ZJD1 = XRADIUS*(abs(ZCLAT0))^(1.-ZRPK)
        ZT1  = (ZJD1)^(2./ZRPK)* (1+ZSLAT0)^2
        ZJD3 = (ZRPK^2*ZRO2)
        ZT2  = ZJD3
        ZT2 = ZT2^(1./ZRPK)
        ZJD1 = (ZT1-ZT2)/(ZT1+ZT2)
        ZJD1 = acos(ZJD1)
        ZJD3 = ZJD1
        zlat = (XPI/2.-ZJD3)/ZRDSDG
        plat(jj,ji)=zlat
      end do
    end do
    
  end if

; Defining the corners of the domain
;====================================
  if (icorner(0,0).eq.icorner(1,1)) then
    icorner(0,0)=0
    icorner(1,0)=JMAX-1
    icorner(0,1)=0
    icorner(1,1)=IMAX-1
  end if
; print ("icorner"+icorner)

  res@mpLimitMode            = "Corners"
  res@mpLeftCornerLatF       = plat(icorner(0,0),icorner(0,1))
  res@mpLeftCornerLonF       = plon(icorner(0,0),icorner(0,1))
  res@mpRightCornerLatF     = plat(icorner(1,0),icorner(1,1))
  res@mpRightCornerLonF     = plon(icorner(1,0),icorner(1,1))

; print ("Corner (0,0); Lat="+res@mpLeftCornerLatF+ \
;                    ", Lon="+res@mpLeftCornerLonF)
; print ("Oppos corner; Lat="+res@mpRightCornerLatF+ \
;                     ", Lon= "+res@mpRightCornerLonF)

;==========================================
; Turn on lat / lon labeling
;==========================================
  res@pmTickMarkDisplayMode = "Always"         ; turn on tickmarks
  res@mpOutlineBoundarySets = "AllBoundaries" ; state boundaries
  res@mpPerimDrawOrder      = "PostDraw"       ; force map perim
;==========================================
; Needed for regional native projection 
;==========================================
  res@tfDoNDCOverlay       = True      
  res@gsnAddCyclic         = False              ; regional data
  
end

;===========================================
;------------------------------------------------------------------------
undef("MESONH_pinter")
function MESONH_pinter( pfield:numeric, loc_param:numeric, ppabs:numeric )
;*************************************************************************
; S. BIELLI
; This is a routine that interpolate  fields on pressure level for plotting
; based on pinter.f90 
; The field to be interpolated must be given at the mass point (grid 1)
; usage : var_inter=MESONHfunction(var_to_interpol, 850., AbsPressure)
; Abs pressure must be in Pa
;

begin

 dimL= dimsizes(loc_param)

; First test for grid = 0

  dimp=dimsizes(ppabs)

  pout=pfield(0:dimL-1,:,:)
  pfield@_FillValue=999
  pout@_FillValue=999
  pout=pout@_FillValue

  do jkp = 0, dimL-1
 	zref=log10(loc_param(jkp)*100.)
	do jloop = 0, dimp(1)-1
  	  do iloop = 0, dimp(2)-1
            kloop=0
            flag=True
	    do while (flag .and. (kloop.lt.(dimp(2)-2)))
	      if (.not.ismissing(ppabs(kloop,jloop,iloop))) then
		zxm=log10(ppabs(kloop,jloop,iloop))
		zxp=log10(ppabs(kloop+1,jloop,iloop))
		if ((zxp-zref)*(zref-zxm) .ge. 0) then
                  pout(jkp,jloop,iloop)= (pfield(kloop,jloop,iloop)*(zxp-zref)+ \
			pfield(kloop+1,jloop,iloop)*(zref-zxm))/ (zxp-zxm)
                  flag=False
                end if
              end if
              kloop=kloop+1
            end do
          end do
        end do
  end do

  return(pout)
	
end

;--------------------------------------------------------------------------------
undef("mnh_map")
function mnh_map(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical)   

begin
;
; This function creates a map plot, and bases the projection on
; the MAP_PROJ attribute in the given file.
;
;   1. Make a copy of the resource list, and set some resources
;      common to all map projections.
;
;   2. Determine the projection being used, and set resources based
;      on that projection.
;
;   3. Create the map plot, and draw and advance the frame
;      (if requested).

    opts = opt_args      ; Make a copy of the resource list
    opts  = True

; Set some resources depending on what kind of map projection is 
; chosen.
;
;   ZRPK != 0 : "Stereographic"
;   ZRPK = 0 : "Mercator"
;=================================================;
; src/mesonh_MOD/mode_gridproj.f90
;=================================================;
  XRADIUS=6371229.0d ; Earth radius (meters)
  XPI=2.0d*asin(1.)    ; Pi
  ZRDSDG= XPI/180.0d          ; Radian to Degree conversion factor
  ZXBM0 = 0.0d
  ZYBM0 = 0.0d

  if(isfilevar(in_file,"RPK"))
    ZRPK=in_file->RPK
    ZLON0=in_file->LON0
    ZLAT0=in_file->LAT0
    ZLATOR=in_file->LATOR
    ZLONOR=in_file->LONOR
    ZBETA=in_file->BETA
  else
   print ("mnh_map: Error no RPK variable in input file")
   return(new(1,graphic))
  end if
 
; Case netcdf from lfi2cdf
   if(isfilevar(in_file,"IMAX"))
       XHAT=in_file->XHAT
       YHAT=in_file->YHAT
       IMAX= dimsizes(XHAT)-2
       JMAX= dimsizes(YHAT)-2
       zdx=XHAT(2)-XHAT(1)
       zdy=YHAT(2)-YHAT(1)
       do ji=0,IMAX-1
         XHAT(ji)=XHAT(ji)+zdx*1.5
       end do
       do jj=0,JMAX-1
         YHAT(jj)=YHAT(jj)+zdy*1.5
       end do
   else
; Case netcdf from extractdia
       XHAT=in_file->W_E_direction
       YHAT=in_file->S_N_direction
       IMAX= dimsizes(XHAT)
       JMAX= dimsizes(YHAT)
   end if
;
  
  lat = new((/JMAX,IMAX/),"double")
  lon = new((/JMAX,IMAX/),"double")


;   Stereographic projection
      if(ZRPK .gt. 0)
        projection          = "Stereographic"
        opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF", 90)
        opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",ZLON0)
        opts@mpCenterRotF  = get_res_value_keep(opts, "mpCenterRotF",ZBETA)
      end if

      if(ZRPK .lt. 0)
        projection          = "Stereographic"
        opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF", -90)
        opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",ZLON0)
        opts@mpCenterRotF  = get_res_value_keep(opts, "mpCenterRotF",ZBETA)
      end if

;   Mercator projection
      if(ZRPK .eq. 0)
        projection          = "Mercator"
      end if

    opts@mpNestTime = get_res_value_keep(opts, "mpNestTime",0)

    
; LAT and LON are not saved in the file 
  if (ZRPK.eq.0) then
    XBETA=0.
    XLAT0=0.   ; map reference latitude  (degrees)
    ZXBM0 = 0.
    ZYBM0 = 0.
    ZCGAM    = cos(-ZRDSDG*XBETA)
    ZSGAM    = sin(-ZRDSDG*XBETA)
    ZRACLAT0 = XRADIUS*cos(ZRDSDG*ZLAT0)
    do ji=0,IMAX-1
      jj=0
      ZXMI0 = XHAT(ji)-ZXBM0
      ZYMI0 = YHAT(jj)-ZYBM0
      zlon  = (ZXMI0*ZCGAM+ZYMI0*ZSGAM)/(ZRACLAT0*ZRDSDG)+ZLONOR
      do jj=0,JMAX-1
        lon(jj,ji)=zlon
      end do
    end do
    do jj=0,JMAX-1
      ji=0
      ZXMI0 = XHAT(ji)-ZXBM0
      ZYMI0 = YHAT(jj)-ZYBM0
      ZT1  = log(tan(XPI/4.+ZLATOR*ZRDSDG/2.))
      ZT2  = (-ZXMI0*ZSGAM+ZYMI0*ZCGAM)/ZRACLAT0
      zlat  = (-XPI/2.+2.*atan(exp(ZT1+ZT2)))/ZRDSDG
      do ji=0,IMAX-1
        lat(jj,ji)=zlat
      end do
    end do
   else
    ZCLAT0  = cos(ZRDSDG*ZLAT0)
    ZSLAT0  = sin(ZRDSDG*ZLAT0)
    ZCLATOR = cos(ZRDSDG*ZLATOR)
    ZSLATOR = sin(ZRDSDG*ZLATOR)
    ZRO0    = (XRADIUS/ZRPK)*(abs(ZCLAT0))^(1.-ZRPK) * \
              ((1.+ZSLAT0)*abs(ZCLATOR)/(1.+ZSLATOR))^ZRPK
    ZGA0    = (ZRPK*(ZLONOR-ZLON0)-ZBETA)*ZRDSDG
    ZXP     = ZXBM0-ZRO0*sin(ZGA0)
    ZYP     = ZYBM0+ZRO0*cos(ZGA0)
    do ji=0,IMAX-1
      do jj=0,JMAX-1
        ZATA = atan2( -(ZXP-XHAT(ji)) , (ZYP-YHAT(jj)) )/ZRDSDG
        zlon  = (ZBETA+ZATA)/ZRPK+ZLON0
        lon(jj,ji)=zlon
        ZRO2 = (XHAT(ji)-ZXP)^2+(YHAT(jj)-ZYP)^2
        ZJD1 = XRADIUS*(abs(ZCLAT0))^(1.-ZRPK)
        ZT1  = (ZJD1)^(2./ZRPK)* (1+ZSLAT0)^2
        ZJD3 = (ZRPK^2*ZRO2)
        ZT2  = ZJD3
        ZT2 = ZT2^(1./ZRPK)
        ZJD1 = (ZT1-ZT2)/(ZT1+ZT2)
        ZJD1 = acos(ZJD1)
        ZJD3 = ZJD1
        zlat = (XPI/2.-ZJD3)/ZRDSDG
        lat(jj,ji)=zlat
      end do
    end do
    end if

      dims = dimsizes(lat)

      do ii = 0, dims(0)-1
      do jj = 0, dims(1)-1
        if ( lon(ii,jj) .lt. 0.0) then
          lon(ii,jj) = lon(ii,jj) + 360.
        end if
      end do
      end do

      opts@start_lat = lat(0,0)
      opts@start_lon = lon(0,0)
      opts@end_lat   = lat(dims(0)-1,dims(1)-1)
      opts@end_lon   = lon(dims(0)-1,dims(1)-1)


; Set some resources common to all map projections.
  opts = set_mp_resources(opts)

      if ( isatt(opts,"ZoomIn") .and. opts@ZoomIn ) then
        y1 = 0
        x1 = 0
        y2 = dims(0)-1
        x2 = dims(1)-1
        if ( isatt(opts,"Ystart") ) then
          y1 = opts@Ystart
          delete(opts@Ystart)
        end if
        if ( isatt(opts,"Xstart") ) then
          x1 = opts@Xstart
          delete(opts@Xstart)
        end if
        if ( isatt(opts,"Yend") ) then
          if ( opts@Yend .le. y2 ) then
            y2 = opts@Yend
          end if
          delete(opts@Yend)
        end if
        if ( isatt(opts,"Xend") ) then
          if ( opts@Xend .le. x2 ) then
            x2 = opts@Xend
          end if
          delete(opts@Xend)
        end if

        opts@mpLeftCornerLatF      = lat(y1,x1)
        opts@mpLeftCornerLonF      = lon(y1,x1)
        opts@mpRightCornerLatF     = lat(y2,x2)
        opts@mpRightCornerLonF     = lon(y2,x2)

        if ( opts@mpRightCornerLonF .lt. 0.0 ) then
          opts@mpRightCornerLonF  = opts@mpRightCornerLonF + 360.0
        end if 

        delete(opts@ZoomIn)
      end if


; The default is not to draw the plot or advance the frame, and
; to maximize the plot in the frame.

  opts@gsnDraw       = get_res_value_keep(opts,"gsnDraw",     False)
  opts@gsnFrame      = get_res_value_keep(opts,"gsnFrame",    False)
  opts@gsnMaximize   = get_res_value_keep(opts,"gsnMaximize", True)

  delete_attrs(opts)                             ; Clean up.
  mp = gsn_map(wks,projection,opts)              ; Create map plot.

  return(mp)                                     ; Return.

end

;--------------------------------------------------------------------------------

undef("mnh_map_overlays")
function mnh_map_overlays(in_file[1]:file, \
                          wks:graphic, \
                          plots[*]:graphic, \
                          opt_arg[1]:logical, \
                          opt_mp[1]:logical) 

; Based on wrf_map_overlays 
; 
; This procedure takes an array of plots and overlays them on a
; base plot - map background.
;
; It will advance the plot and cleanup, unless you set the
; PanelPlot resource to True.
;
;  Attributes recognized by this procedure:
;     FramePlot
;     PanelPlot
;     NoTitles                  (don't do any titles) 
;     CommonTitle & PlotTile is used to overwrite field titles
;        CommonTitle will super-seed NoTitles
;
; If FramePlot False, then Draw the plot but do not Frame.
; In this case a user want to add to the drawing, and will
; have to advance the Frame manually in the script.
;
; If the "NoTitles" attribute exists and is set True, then
; don't create the top-left titles, and leave the main titles alone.
; This resource can be useful if you are planning to panel
; the plots.
;
; If PanelPlot is set to True, then this flags to wrf_map_overlays
; that these plots are going to be eventually paneled (likely
; by gsn_panel), and hence 1) draw and frame should not be called
; (unless gsnDraw and/or gsnFrame are explicitly set to True),
; and 2) the overlays and titles should not be removed with
; NhlRemoveOverlay and NhlRemoveAnnotation.
;
begin

  opts     = opt_arg         ; Make a copy of the resource lists
  opt_mp_2 = opt_mp      

  ; Let's make the map first
  base = mnh_map(wks,in_file,opt_mp_2)

  no_titles  = get_res_value(opts,"NoTitles",False)     ; Do we want field titles?
  com_title  = get_res_value(opts,"CommonTitle",False)     ; Do we have a common title?
  if ( com_title ) then
    plot_title = get_res_value(opts,"PlotTitle","  ")
    no_titles = True
  end if
  
  call_draw  = True
  call_frame = get_res_value(opts,"FramePlot",True)     ; Do we want to frame the plot?
  panel_plot = get_res_value(opts,"PanelPlot",False)    ; Are we paneling?
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  nplots = dimsizes(plots)
;  font_color = "Black"

  do i=0,nplots-1
    if(.not.ismissing(plots(i))) then
;      class_name = NhlClassName(plots(i))
;      print(class_name)
;      if(class_name.eq."contourPlotClass") then
;        getvalues plots(i)
;          "cnFillOn"    : fill_on
;          "cnLineColor" : line_color
;        end getvalues
;        if (.not.fill_on) then
;          font_color = line_color      
;        end if       
;      end if
      if(.not.no_titles) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        if(i.eq.0) then
          SubTitles = SubTitle
        else
          SubTitles = SubTitles + "~C~" + SubTitle
        end if
      end if
      if(com_title .and. i .eq. nplots-1) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        SubTitles = plot_title
      end if
      setvalues plots(i)
        "tfDoNDCOverlay" : True
        "tiMainOn"       : False
      end setvalues
      overlay(base,plots(i))
    else
      print("mnh_map_overlays: Warning: overlay plot #" + i + " is not valid.")
    end if
  end do

  if(.not.no_titles .or. com_title) then
    font_height = get_res_value_keep(opts,"FontHeightF",0.01)
    txt = create "map_titles" textItemClass wks
      "txString"      : SubTitles
      "txFontHeightF" : font_height
     ;"txFontColor"   : font_color
    end create
    anno = NhlAddAnnotation(base,txt)
    setvalues anno
      "amZone"           : 3
      "amJust"           : "BottomLeft"
      "amSide"           : "Top"
      "amParallelPosF"   : 0.005
      "amOrthogonalPosF" : 0.03
      "amResizeNotify"   : False
    end setvalues
    base@map_titles = anno
  end if
;
; gsnDraw and gsnFrame default to False if panel plot.
;
  if(panel_plot) then
    call_draw = False
    call_frame= False
  end if


  opts@gsnDraw     = get_res_value_keep(opts,"gsnDraw",     call_draw)
  opts@gsnFrame    = get_res_value_keep(opts,"gsnFrame",    call_frame)

  draw_and_frame(wks,base,opts@gsnDraw,opts@gsnFrame,False, \
                 opts@gsnMaximize)

  if(.not.panel_plot) then
    do i=0,nplots-1
      if(.not.ismissing(plots(i))) then
        NhlRemoveOverlay(base,plots(i),False)
      else
        print("wrf_remove_map_overlays: Warning: overlay plot #" + i + " is not valid.")
        print("                        Nothing to remove.")
      end if
    end do
  end if

  if(.not.no_titles.and..not.panel_plot) then
    if(isatt(base,"map_titles")) then
      NhlRemoveAnnotation(base,base@map_titles)
      delete(base@map_titles)
    end if
  end if

return(base)
end

;--------------------------------------------------------------------------------
undef("wrf_user_intrp3d")
function wrf_user_intrp3d( var3d:numeric, z_in:numeric, \
                           plot_type:string, \
                           loc_param:numeric, angle:numeric, opts:logical )

; var3d      - 3d field to interpolate (all input fields must be unstaggered)
; z_in       - interpolate to this field (either p/z)
; plot_type  - interpolate horizontally "h", or vertically "v"
; loc_param  - level(s) for horizontal plots (eg. 500hPa ; 3000m - scalar),
;              plane for vertical plots (2 values representing an xy point
;              on the model domain through which the vertical plane will pass
;              OR 4 values specifying start and end values
; angle      - 0.0 for horizontal plots, and
;              an angle for vertical plots - 90 represent a WE cross section
; opts         Used IF opts is TRUE, else use loc_param and angle to determine crosssection

begin


  if(plot_type .eq. "h" ) then   ;  horizontal cross section needed

     dimL = dimsizes(loc_param)

     dims = dimsizes(var3d)
     nd = dimsizes(dims)

     dimX = dims(nd-1)
     dimY = dims(nd-2)
     dimZ = dims(nd-3)
     dim4 = 1
     dim5 = 1
     if ( nd .eq. 4 ) then
       dim4 = dims(nd-4)
     end if
     if ( nd .eq. 5 ) then
       dim4 = dims(nd-4)
       dim5 = dims(nd-5)
     end if

     var3  = new ( (/ dim5, dim4, dimZ, dimY, dimX /) , typeof(var3d) )
     z     = new ( (/ dim5, dim4, dimZ, dimY, dimX /) , typeof(var3d) )
     var2d = new ( (/ dim5, dim4, dimL, dimY, dimX /) , typeof(var3d) )

     if ( nd .eq. 5 ) then
       var3 = var3d
          z =  z_in
     end if
     if ( nd .eq. 4 ) then
       var3(0,:,:,:,:) = var3d(:,:,:,:)
          z(0,:,:,:,:) =  z_in(:,:,:,:)
     end if
     if ( nd .eq. 3 ) then
       var3(0,0,:,:,:) = var3d(:,:,:)
          z(0,0,:,:,:) =  z_in(:,:,:)
     end if


     if ( z(0,0,0,0,0) .gt. 500.) then
       ; We must be interpolating to pressure
       ; This routine needs input field and level in hPa - lets make sure of this
       if ( z(0,0,0,0,0) .gt. 2000. ) then
         ; looks like we have Pa as input - make this hPa
         z = z * 0.01
       end if
       if ( loc_param(0) .gt. 2000. ) then
         ; looks like the input was specified in Pa - change this
         loc_param = loc_param * 0.01
       end if
     end if

     do il = 0,dimL-1
       var = wrf_interp_3d_z(var3,z,loc_param(il))
       var2d(:,:,il,:,:) = var(:,:,:,:)
     end do

     copy_VarAtts(var3d,var3)
     if(isatt(var3,"description")) then
       delete_VarAtts(var3,(/"description"/))
     end if
     if(isatt(var3,"units")) then
       delete_VarAtts(var3,(/"units"/))
     end if
     if(isatt(var3,"MemoryOrder")) then
       delete_VarAtts(var3,(/"MemoryOrder"/))
     end if
     if(isatt(var3,"_FillValue")) then
       delete_VarAtts(var3,(/"_FillValue"/))
     end if
     copy_VarAtts(var3,var2d)

     nn = nd-2
     var2d!nn = "plevs"

     if ( dimL .gt. 1 ) then
       if ( nd .eq. 5 ) then
         return( var2d )
       end if
       if ( nd .eq. 4 ) then
         return( var2d(0,:,:,:,:) )
       end if
       if ( nd .eq. 3 ) then
         return( var2d(0,0,:,:,:) )
       end if
     else
       if ( z(0,0,0,0,0) .gt. 500.) then
          var2d@PlotLevelID = loc_param + " hPa"
       else
          var2d@PlotLevelID = .001*loc_param + " km"
       end if
       if ( nd .eq. 5 ) then
         return( var2d(:,:,0,:,:) )
       end if
       if ( nd .eq. 4 ) then
         return( var2d(0,:,0,:,:) )
       end if
       if ( nd .eq. 3 ) then
         return( var2d(0,0,0,:,:) )
       end if
     end if


  end if




  if(plot_type .eq. "v" ) then   ;  vertical cross section needed

     dims = dimsizes(var3d)
     if ( dimsizes(dims) .eq. 4 ) then
       if ( z_in(0,0,0,0) .gt. 500.) then
         ; We must be interpolating to pressure
           ; This routine needs input field and level in hPa - lets make sure of this
         if ( z_in(0,0,0,0) .gt. 2000. ) then
           ; looks like we have Pa as input - make this hPa
           z_in = z_in * 0.01
         end if
       end if
       z = z_in(0,:,:,:)
     else
       if ( z_in(0,0,0) .gt. 500.) then
         ; We must be interpolating to pressure
           ; This routine needs input field and level in hPa - lets make sure of this
         if ( z_in(0,0,0) .gt. 2000. ) then
           ; looks like we have Pa as input - make this hPa
           z_in = z_in * 0.01
         end if
       end if
       z = z_in
     end if

; set vertical cross section
     if (opts) then
       xy = wrf_user_set_xy( z, loc_param(0)-1, loc_param(1)-1, \    ; the -1 is for NCL dimensions
                                loc_param(2)-1, loc_param(3)-1, \
                                angle, opts )
     else
       xy = wrf_user_set_xy( z, loc_param(0), loc_param(1), \
                                0.0, 0.0, angle, opts )
     end if
     xp = dimsizes(xy)


; first we interp z
     var2dz   = wrf_interp_2d_xy( z, xy)

;  interp to constant z grid
     if(var2dz(0,0) .gt. var2dz(1,0) ) then  ; monotonically decreasing coordinate
        z_max = floor(max(z)/10)*10     ; bottom value
        z_min = ceil(min(z)/10)*10      ; top value
        dz = 1.
        nlevels = tointeger( (z_max-z_min)/dz)
        z_var2d = new( (/nlevels/), typeof(z))
        z_var2d(0) = z_max
        dz = -dz
     else
        z_max = max(z)
        z_min = 0.
;; MODI SOLINE
;        dz = 0.01 * z_max
        dz = 0.001 * z_max
        nlevels = tointeger( z_max/dz )
        z_var2d = new( (/nlevels/), typeof(z))
        z_var2d(0) = z_min
     end if
;     print("nlevels="+nlevels)
;     print("dz="+dz)

     do i=1, nlevels-1
        z_var2d(i) = z_var2d(0)+i*dz
     end do


; interp the variable
     if ( dimsizes(dims) .eq. 4 ) then
       var2d = new( (/dims(0), nlevels, xp(0)/), typeof(var2dz))
       do it = 0,dims(0)-1
         var2dtmp = wrf_interp_2d_xy( var3d(it,:,:,:), xy)
         do i=0,xp(0)-1
            var2d(it,:,i) = wrf_interp_1d( var2dtmp(:,i), var2dz(:,i), z_var2d)
         end do
       end do
       var2d!0 = var3d!0
       var2d!1 = "Vertical"
       var2d!2 = "Horizontal"
     else
       var2d = new( (/nlevels, xp(0)/), typeof(var2dz))
       var2dtmp = wrf_interp_2d_xy( var3d, xy)
       do i=0,xp(0)-1
          var2d(:,i) = wrf_interp_1d( var2dtmp(:,i), var2dz(:,i), z_var2d)
       end do
       var2d!0 = "Vertical"
       var2d!1 = "Horizontal"
     end if


     st_x = tointeger(xy(0,0)) + 1
     st_y = tointeger(xy(0,1)) + 1
     ed_x = tointeger(xy(xp(0)-1,0)) + 1
     ed_y = tointeger(xy(xp(0)-1,1)) + 1
     if (opts) then
       var2d@Orientation = "Cross-Sesion: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ")"
     else
       var2d@Orientation = "Cross-Sesion: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ") ; center=(" + \
                            loc_param(0) + "," + loc_param(1) + \
                            ") ; angle=" + angle
     end if

     return(var2d)
end if


end

