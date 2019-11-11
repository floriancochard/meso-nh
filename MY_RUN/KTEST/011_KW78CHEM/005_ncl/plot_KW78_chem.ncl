;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl" 
; ================================================;
;================================================;
begin
;=================================================;
; open file and read in data
;=================================================;
  fichier1 = addfile("KWRAI.1.SEGCH.004.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_KW78_chem")
;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = fichier1->JPHEXT
  mdims = getfilevardimsizes(fichier1,"UT") ; get dimension sizes
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2
;=================================================;
; lecture des différents champs
;=================================================;
zs  = fichier1->ZS(jphext:jmax+jphext-1,jphext:imax+jphext-1) ; ZS
zhat  = fichier1->ZHAT(1:kmax+1) ; ZHAT
yhat  = fichier1->YHAT(1:jmax+jphext-1+1) ; YHAT
xhat  = fichier1->XHAT(1:imax+jphext-1+1) ; XHAT

wt= fichier1->WT(0,1:kmax+1,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; WT
vt= fichier1->VT(0,1:kmax,jphext:jmax+jphext-1+1,jphext:imax+jphext-1) ; VT
ut= fichier1->UT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1+1) ; UT
tht1= fichier1->THT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; THT
tht1@long_name="Potential Temperature"
tht1@units="K"
  lsthm1= fichier1->LSTHM(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
lsthm1@long_name="Large SCale Potential Temperature"
lsthm1@units="K"

pabst1= fichier1->PABST(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; PABST
pabst1@long_name="Pressure"
pabst1@units="Pa"
rvt1= fichier1->RVT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RVT
rvt1@long_name="Vapor mixing ratio"
rvt1@units="kg/kg"
rct1= fichier1->RCT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RCT
rct1@long_name="Cloud mixing ratio"
rct1@units="kg/kg"
rrt1= fichier1->RRT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RRT
rrt1@long_name="Rain mixing ratio"
rrt1@units="kg/kg"
rst1= fichier1->RST(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RRT
rst1@long_name="Snow mixing ratio"
rst1@units="kg/kg"
rgt1= fichier1->RGT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; RRT
rgt1@long_name="Graupel mixing ratio"
rgt1@units="kg/kg"

inprr= fichier1->INPRR(0,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
inprr@long_name="Instantaneous Precipitation Rain Rate"
inprr@units="mm/h"
acprr= fichier1->ACPRR(0,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
acprr@long_name="Accumulated Precipitation Rain Rate"
acprr@units="mm"


cot1= fichier1->COT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; COT
cot1@long_name="carbon monoxide"
cot1@units="ppp"
o3t1= fichier1->O3T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; O3T
o3t1@long_name="Ozone"
o3t1@units="ppp"

hchot1= fichier1->HCHOT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; HCHOT
hchot1@long_name="HCHO"
hchot1@units="ppp"
wc_hchot1= fichier1->WC_HCHOT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
wc_hchot1@long_name="WC_HCHO"
wc_hchot1@units="ppp"
wr_hchot1= fichier1->WR_HCHOT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
wr_hchot1@long_name="WR_HCHO"
wr_hchot1@units="ppp"
ic_hchot1= fichier1->IC_HCHOT(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
ic_hchot1@long_name="IC_HCHO"
ic_hchot1@units="ppp"

so2t1= fichier1->SO2T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; SO2T
so2t1@long_name="SO2"
so2t1@units="ppp"
wc_so2t1= fichier1->WC_SO2T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
wc_so2t1@long_name="WC_SO2"
wc_so2t1@units="ppp"
wr_so2t1= fichier1->WR_SO2T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
wr_so2t1@long_name="WR_SO2"
wr_so2t1@units="ppp"
ic_so2t1= fichier1->IC_SO2T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
ic_so2t1@long_name="IC_SO2"
ic_so2t1@units="ppp"

hno3t1= fichier1->HNO3T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; HNO3T
hno3t1@long_name="HNO3"
hno3t1@units="ppp"
wc_hno3t1= fichier1->WC_HNO3T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
wc_hno3t1@long_name="WC_HNO3"
wc_hno3t1@units="ppp"
wr_hno3t1= fichier1->WR_HNO3T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
wr_hno3t1@long_name="WR_HNO3"
wr_hno3t1@units="ppp"
ic_hno3t1= fichier1->IC_HNO3T(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)
ic_hno3t1@long_name="IC_HNO3"
ic_hno3t1@units="ppp"

phc1= fichier1->PHC(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) 
phc1@long_name="PHC"
phc1@units=" "

phr1= fichier1->PHR(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1)  
phr1@long_name="PHR"
phr1@units=" "

;YYYYDDMM1= fichier1->DTCUR__TDATE
SS1=fichier1->DTCUR

;=================================================;
; Récupération de la date 
;=================================================;
;year1=YYYYDDMM1(0)
;day1=YYYYDDMM1(1)
;min1=YYYYDDMM1(2)
time1=SS1
;=================================================;
; On mets toutes les variables sur la grille 1 
;=================================================;
; rq : cas cartésien
ut1 = wrf_user_unstagger(ut,"X")
ut1@long_name="Zonal wind"
ut1@units="m/s"
vt1 = wrf_user_unstagger(vt,"Y")
vt1@long_name="Meridian wind"
vt1@units="m/s"
wt1 = wrf_user_unstagger(wt,"Z")
wt1@long_name="Vertical wind"
wt1@units="m/s"

;=================================================;
; Set some other basic resources
;=================================================;
  resmap = True
  resmap@gsnFrame = False
  resmap@gsnDraw = False
  resmap@gsnMaximize = True
  resmap@gsnPaperOrientation = "portrait" 
  resmap@gsnSpreadColors    	= True       	; use full range of colormap
  resmap@tiYAxisString =" "
  resmap@cnFillOn               = True ; turn on color fill
  resmap@cnLinesOn               = False ; turn off contour lines
;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
    nzh=new(kmax,double)
    do k=0,kmax-1
     nzh(k)=(zhat(k)+zhat(k+1))/2.
    end do

; Create a (altitude des niveaux modèle)
    alt=new(dimsizes(tht1),double)
    zcoef=new((/jmax,imax/),double)
    zcoef=1.-zs/nzh(kmax-1)

    do i=0,imax-1
      do j=0,jmax-1
       alt(:,j,i) = nzh*zcoef(j,i)+zs(j,i)
      end do
    end do


;=================================================;
; TRACE
;=================================================;
  gsn_define_colormap(wks,"rainbow") ; Choose colormap

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 0.25
opts@cnMinLevelValF    = 0
opts@cnMaxLevelValF    = 6

plot_inprr = gsn_csm_contour(wks,inprr,opts)
draw(plot_inprr)
frame(wks)

opts@cnLevelSpacingF    = 0.005
opts@cnMinLevelValF    = 0
plot_acprr = gsn_csm_contour(wks,acprr,opts)
draw(plot_acprr)
frame(wks)

res=resmap
res@gsnScalarContour=True
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF    = 5
res@cnMinLevelValF    = 98250
res@cnMaxLevelValF    = 98350
res@vcMinDistanceF = 0.03 ; espacement des fleches
res@vcRefLengthF=0.1      ; taille des flèches
plot=gsn_csm_vector_scalar(wks,ut1(0,:,:),vt1(0,:,:),pabst1(0,:,:),res)
draw(plot)
frame(wks)


rct_plane = wrf_user_intrp3d(rct1,alt,"h",3000,0.,False)
ut_plane = wrf_user_intrp3d(ut1,alt,"h",3000,0.,False)
vt_plane = wrf_user_intrp3d(vt1,alt,"h",3000,0.,False)

res2=resmap
res2@gsnScalarContour=True
res2@vcMinDistanceF = 0.04 ; espacement des fleches
res2@vcRefLengthF=0.1      ; taille des flèches
res2@cnLevelSpacingF    = 0.00005
res2@cnMinLevelValF    = 0
res2@tiXAxisString = "Z=3000m"
plot=gsn_csm_vector_scalar(wks,ut_plane(:,:),vt_plane(:,:),rct_plane(:,:),res2)
draw(plot)
frame(wks)

rct5000_plane = wrf_user_intrp3d(rct1,alt,"h",5000,0.,False)
ut5000_plane = wrf_user_intrp3d(ut1,alt,"h",5000,0.,False)
vt5000_plane = wrf_user_intrp3d(vt1,alt,"h",5000,0.,False)
res2@tiXAxisString = "Z=5000m"
plot=gsn_csm_vector_scalar(wks,ut5000_plane(:,:),vt5000_plane(:,:),rct5000_plane(:,:),res2)
draw(plot)
frame(wks)


cot3000_plane = wrf_user_intrp3d(cot1,alt,"h",3000,0.,False)
o3t3000_plane = wrf_user_intrp3d(o3t1,alt,"h",3000,0.,False)
cot5000_plane = wrf_user_intrp3d(cot1,alt,"h",5000,0.,False)
o3t5000_plane = wrf_user_intrp3d(o3t1,alt,"h",5000,0.,False)
res@cnLevelSelectionMode = "AutomaticLevels"
res2@tiXAxisString = "Z=3000m"
plot=gsn_csm_contour(wks,cot3000_plane(:,:),res2)
draw(plot)
frame(wks)
plot=gsn_csm_contour(wks,o3t3000_plane(:,:),res2)
draw(plot)
frame(wks)
res2@tiXAxisString = "Z=5000m"
plot=gsn_csm_contour(wks,cot5000_plane(:,:),res2)
draw(plot)
frame(wks)
plot=gsn_csm_contour(wks,o3t5000_plane(:,:),res2)
draw(plot)
frame(wks)


;============================================
;   COUPE VERTICALE
;============================================
; definition de la coupe verticale
opt_cv=False
plane = new(2,float)
plane =(/3,0/) ; point de depart attention i puis j
angle=43 ; attention 0=S-N et 90=O-E


rvt_cv=wrf_user_intrp3d(rvt1,alt,"v",plane,angle,False)
rct_cv=wrf_user_intrp3d(rct1,alt,"v",plane,angle,False)
rrt_cv=wrf_user_intrp3d(rrt1,alt,"v",plane,angle,False)
rst_cv=wrf_user_intrp3d(rst1,alt,"v",plane,angle,False)
rgt_cv=wrf_user_intrp3d(rgt1,alt,"v",plane,angle,False)

ut_cv=wrf_user_intrp3d(ut1,alt,"v",plane,angle,False)
vt_cv=wrf_user_intrp3d(vt1,alt,"v",plane,angle,False)
wt_cv=wrf_user_intrp3d(wt1,alt,"v",plane,angle,False)


rvt_cv@long_name="Vapor mixing ratio"
rct_cv@long_name="Cloud mixing ratio"
rrt_cv@long_name="Rain mixing ratio"
rgt_cv@long_name="Graupel mixing ratio"
rst_cv@long_name="Snow mixing ratio"
rgst_cv=rgt_cv+rst_cv
rgst_cv@long_name="Graupel + Snow mixing ratio"


; projection du vent horizontal dans la coupe verticale
wind_cv=ut_cv*sin(angle*3.14/180.)+vt_cv*cos(angle*3.14/180)


; projection chimie

cot_cv=wrf_user_intrp3d(cot1,alt,"v",plane,angle,False)
hchot_cv=wrf_user_intrp3d(hchot1,alt,"v",plane,angle,False)
wc_hchot_cv=wrf_user_intrp3d(wc_hchot1,alt,"v",plane,angle,False)
wr_hchot_cv=wrf_user_intrp3d(wr_hchot1,alt,"v",plane,angle,False)
ic_hchot_cv=wrf_user_intrp3d(ic_hchot1,alt,"v",plane,angle,False)
so2t_cv=wrf_user_intrp3d(so2t1,alt,"v",plane,angle,False)
wc_so2t_cv=wrf_user_intrp3d(wc_so2t1,alt,"v",plane,angle,False)
wr_so2t_cv=wrf_user_intrp3d(wr_so2t1,alt,"v",plane,angle,False)
ic_so2t_cv=wrf_user_intrp3d(ic_so2t1,alt,"v",plane,angle,False)
hno3t_cv=wrf_user_intrp3d(hno3t1,alt,"v",plane,angle,False)
wc_hno3t_cv=wrf_user_intrp3d(wc_hno3t1,alt,"v",plane,angle,False)
wr_hno3t_cv=wrf_user_intrp3d(wr_hno3t1,alt,"v",plane,angle,False)
ic_hno3t_cv=wrf_user_intrp3d(ic_hno3t1,alt,"v",plane,angle,False)
phc_cv=wrf_user_intrp3d(phc1,alt,"v",plane,angle,False)
phr_cv=wrf_user_intrp3d(phr1,alt,"v",plane,angle,False)


opt_cv=resmap
opt_cv@gsnScalarContour=True
opt_cv@vcMinDistanceF = 0.03 ; espacement des fleches
opt_cv@vcRefLengthF=0.1      ; taille des flèches
opt_cv@cnLevelSpacingF    = 0.0005
opt_cv@cnMinLevelValF    = 0
plot=gsn_csm_vector_scalar(wks,wind_cv(:,:),wt_cv(:,:),rvt_cv(:,:),opt_cv)
draw(plot)
frame(wks)


res_rc=resmap
res_rc@cnFillOn               = True ; turn on color fill
res_rc@cnLinesOn               = False ; turn off contour lines
res_rc@cnLevelSpacingF    = 0.00005
res_rc@cnMinLevelValF    = 0.
res_rc@gsnLeftString =" "
res_rc@tiXAxisString="Cloud mixing ratio (kg/kg)"
plot_rc=gsn_csm_contour(wks,rct_cv,res_rc)

res_rr=resmap
res_rr@cnFillOn               = False ; turn off color fill
res_rr@cnLinesOn               = True ; turn on contour lines
res_rr@cnLevelSpacingF    = 0.00025
res_rr@cnMinLevelValF    = 0.00025
plot_rr=gsn_csm_contour(wks,rrt_cv,res_rr)

overlay(plot_rc,plot_rr)
draw(plot_rc)
frame(wks)

tht_cv=wrf_user_intrp3d(tht1,alt,"v",plane,angle,False)
lsthm_cv=wrf_user_intrp3d(lsthm1,alt,"v",plane,angle,False)
res_th=resmap
res_th@cnLevelSpacingF    = 0.2
res_th@cnMinLevelValF    = -4
res_th@cnMaxLevelValF    = 4
res_th@tiMainString="THT-LSTHM (K)"
gsn_define_colormap(wks,"BlRe") ; Choose colormap
plot_th=gsn_csm_contour(wks,tht_cv-lsthm_cv,res_th)
draw(plot_th)
frame(wks)

res_wt=resmap
res_wt@cnLevelSpacingF    = 0.5
res_wt@cnMinLevelValF    = -4
res_wt@cnMaxLevelValF    = 10
res_wt@tiMainString="WT"
plot_wt=gsn_csm_contour(wks,wt_cv,res_wt)
draw(plot_wt)
frame(wks)


  gsn_define_colormap(wks,"rainbow") ; Choose colormap

res_cot=resmap
res_cot@tiMainString="Carbon Monoxyde"
plot_cot=gsn_csm_contour(wks,cot_cv,res_cot)

res_hchot=resmap
res_hchot@tiMainString="HCHO"
plot_hchot=gsn_csm_contour(wks,hchot_cv,res_hchot)

res_wc_hchot=resmap
res_wc_hchot@gsnLeftString =" "
res_wc_hchot@gsnRightString =" "
res_wc_hchot@tiMainString="WC_HCHO (ppp) + Cloud mixing ratio (kg/kg) "
plot_wc_hchot=gsn_csm_contour(wks,wc_hchot_cv,res_wc_hchot)

res_wr_hchot=resmap
res_wr_hchot@gsnLeftString =" "
res_wr_hchot@gsnRightString =" "
res_wr_hchot@tiMainString="WR_HCHO (ppp) + rain mixing ratio (kg/kg) "
plot_wr_hchot=gsn_csm_contour(wks,wr_hchot_cv,res_wr_hchot)
res_ic_hchot=resmap
res_ic_hchot@gsnLeftString =" "
res_ic_hchot@gsnRightString =" "
res_ic_hchot@tiMainString="IC_HCHO (ppp) + snow and graupel mixing ratio (kg/kg) "
plot_ic_hchot=gsn_csm_contour(wks,ic_hchot_cv,res_ic_hchot)



res_rc2=resmap
res_rc2@cnFillOn               = False ; turn on color fill
res_rc2@cnLinesOn               = True ; turn off contour lines
res_rc2@cnLevelSpacingF    = 0.00005
res_rc2@cnMinLevelValF    = 0.
res_rc2@gsnLeftString =" "
res_rc2@gsnRightString =" "
plot_rc2=gsn_csm_contour(wks,rct_cv,res_rc2)


res_rr2=resmap
res_rr2@cnFillOn               = False ; turn off color fill
res_rr2@cnLinesOn               = True ; turn on contour lines
res_rr2@cnLevelSpacingF    = 0.00025
res_rr2@cnMinLevelValF    = 0.00025
res_rr2@gsnLeftString =" "
res_rr2@gsnRightString =" "
plot_rr2=gsn_csm_contour(wks,rrt_cv,res_rr2)

res_rgs2=resmap
res_rgs2@cnFillOn               = False ; turn on color fill
res_rgs2@cnLinesOn               = True ; turn off contour lines
res_rgs2@gsnLeftString =" "
res_rgs2@gsnRightString =" "
plot_rgs2=gsn_csm_contour(wks,rgst_cv,res_rgs2)


draw(plot_hchot)
frame(wks)
overlay(plot_wc_hchot,plot_rc2)
draw(plot_wc_hchot)
frame(wks)
overlay(plot_wr_hchot,plot_rr2)
draw(plot_wr_hchot)
frame(wks)
overlay(plot_ic_hchot,plot_rgs2)
draw(plot_ic_hchot)
frame(wks)

res_so2t=resmap
res_so2t@tiMainString="SO2"
plot_so2t=gsn_csm_contour(wks,so2t_cv,res_so2t)

res_wc_so2t=resmap
res_wc_so2t@gsnLeftString =" "
res_wc_so2t@gsnRightString =" "
res_wc_so2t@tiMainString="WC_SO2 (ppp) + Cloud mixing ratio (kg/kg) "
plot_wc_so2t=gsn_csm_contour(wks,wc_so2t_cv,res_wc_so2t)

res_wr_so2t=resmap
res_wr_so2t@gsnLeftString =" "
res_wr_so2t@gsnRightString =" "
res_wr_so2t@tiMainString="WR_SO2 (ppp) + rain mixing ratio (kg/kg) "
plot_wr_so2t=gsn_csm_contour(wks,wr_so2t_cv,res_wr_so2t)
res_ic_so2t=resmap
res_ic_so2t@gsnLeftString =" "
res_ic_so2t@gsnRightString =" "
res_ic_so2t@tiMainString="IC_SO2 (ppp) + snow and graupel mixing ratio (kg/kg) "
plot_ic_so2t=gsn_csm_contour(wks,ic_so2t_cv,res_ic_so2t)

plot_rc3=gsn_csm_contour(wks,rct_cv,res_rc2)
plot_rr3=gsn_csm_contour(wks,rrt_cv,res_rr2)
plot_rgs3=gsn_csm_contour(wks,rgst_cv,res_rgs2)


draw(plot_so2t)
frame(wks)
overlay(plot_wc_so2t,plot_rc3)
draw(plot_wc_so2t)
frame(wks)
overlay(plot_wr_so2t,plot_rr3)
draw(plot_wr_so2t)
frame(wks)
overlay(plot_ic_so2t,plot_rgs3)
draw(plot_ic_so2t)
frame(wks)

res_hno3t=resmap
res_hno3t@tiMainString="HNO3"
plot_hno3t=gsn_csm_contour(wks,hno3t_cv,res_hno3t)

res_wc_hno3t=resmap
res_wc_hno3t@gsnLeftString =" "
res_wc_hno3t@gsnRightString =" "
res_wc_hno3t@tiMainString="WC_HNO3 (ppp) + Cloud mixing ratio (kg/kg) "
plot_wc_hno3t=gsn_csm_contour(wks,wc_hno3t_cv,res_wc_hno3t)
res_wr_hno3t=resmap
res_wr_hno3t@gsnLeftString =" "
res_wr_hno3t@gsnRightString =" "
res_wr_hno3t@tiMainString="WR_HNO3 (ppp) + rain mixing ratio (kg/kg) "
plot_wr_hno3t=gsn_csm_contour(wks,wr_hno3t_cv,res_wr_hno3t)
res_ic_hno3t=resmap
res_ic_hno3t@gsnLeftString =" "
res_ic_hno3t@gsnRightString =" "
res_ic_hno3t@tiMainString="IC_HNO3 (ppp) + snow and graupel mixing ratio (kg/kg) "
plot_ic_hno3t=gsn_csm_contour(wks,ic_hno3t_cv,res_ic_hno3t)

plot_rc4=gsn_csm_contour(wks,rct_cv,res_rc2)
plot_rr4=gsn_csm_contour(wks,rrt_cv,res_rr2)
plot_rgs4=gsn_csm_contour(wks,rgst_cv,res_rgs2)


draw(plot_hno3t)
frame(wks)
overlay(plot_wc_hno3t,plot_rc4)
draw(plot_wc_hno3t)
frame(wks)
overlay(plot_wr_hno3t,plot_rr4)
draw(plot_wr_hno3t)
frame(wks)
overlay(plot_ic_hno3t,plot_rgs4)
draw(plot_ic_hno3t)
frame(wks)

res_phc=resmap
res_phc@tiMainString="PHC"
plot_phc=gsn_csm_contour(wks,phc_cv,res_phc)
draw(plot_phc)
frame(wks)

res_phr=resmap
res_phr@tiMainString="PHR"
plot_phr=gsn_csm_contour(wks,phr_cv,res_phr)
draw(plot_phr)
frame(wks)


end
