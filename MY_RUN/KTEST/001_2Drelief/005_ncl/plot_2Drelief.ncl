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
  fichier1 = addfile("EXPER.1.HYD2D.002.nc", "r")
  fichier2 = addfile("EXPER.1.HYD2D.003.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_2Drelief")
;=================================================;
; lecture des différents champs
;=================================================;
  jphext  = fichier1->JPHEXT
  mdims = getfilevardimsizes(fichier1,"UT") 
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2


zs  = fichier1->ZS(0,jphext:imax+jphext-1) ; ZS
zhat  = fichier1->ZHAT(1:kmax+1) ; ZHAT
xhat  = fichier1->XHAT(jphext:imax+jphext-1+1) ; XHAT

lsthm1 = fichier1->LSTHM(0,1:kmax,0,jphext:imax+jphext-1) ; LSTHM
lsthm1@long_name="LSTHM"
lsthm1@units="K"

lsum1_old = fichier1->LSUM(0,1:kmax,0,jphext:imax+jphext-1+1) ; LSUM
wt1_old= fichier1->WT(0,1:kmax+1,0,jphext:imax+jphext-1) ; WT
ut1_old= fichier1->UT(0,1:kmax,0,jphext:imax+jphext-1+1) ; UT

tht1= fichier1->THT(0,1:kmax,0,jphext:imax+jphext-1) ; THT
tht1@long_name="Potential Temperature"
tht1@units="K"

rvt1= fichier1->RVT(0,1:kmax,0,jphext:imax+jphext-1) ; RVT
rvt1@long_name="Vapor mixing ratio"
rvt1@units="g/kg"

cflu1= fichier1->CFLU(0,1:kmax,0,jphext:imax+jphext-1) ; CFLU
cflu1@long_name="CFLU"
cflu1@units="-"

cflw1= fichier1->CFLW(0,1:kmax,0,jphext:imax+jphext-1) ; CFLW
cflw1@long_name="CFLW"
cflw1@units=" "

;YYYYDDMM1= fichier1->DTCUR__TDATE
SS1=fichier1->DTCUR
;YYYYDDMM2= fichier2->DTCUR__TDATE
SS2=fichier2->DTCUR

lsthm2 = fichier2->LSTHM(0,1:kmax,0,jphext:imax+jphext-1) ; LSTHM
lsthm2@long_name="LSTHM"
lsthm2@units="K"

lsum2_old = fichier2->LSUM(0,1:kmax,0,jphext:imax+jphext-1+1) ; LSUM
wt2_old= fichier2->WT(0,1:kmax+1,0,jphext:imax+jphext-1) ; WT
ut2_old= fichier2->UT(0,1:kmax,0,jphext:imax+jphext-1+1) ; UT

tht2= fichier2->THT(0,1:kmax,0,jphext:imax+jphext-1) ; THT
tht2@long_name="Potential Temperature"
tht2@units="K"

rvt2= fichier2->RVT(0,1:kmax,0,jphext:imax+jphext-1) ; RVT
rvt2@long_name="Vapor mixing ratio"
rvt2@units="g/kg"

cflu2= fichier2->CFLU(0,1:kmax,0,jphext:imax+jphext-1) ; CFLU
cflu2@long_name="CFLU"
cflu2@units="-"

cflw2= fichier2->CFLW(0,1:kmax,0,jphext:imax+jphext-1) ; CFLW
cflw2@long_name="CFLW"
cflw2@units=" "

;=================================================;
; Récupération de la date 
;=================================================;
;year1=YYYYDDMM1(0)
;day1=YYYYDDMM1(1)
;min1=YYYYDDMM1(2)
time1=SS1
;year2=YYYYDDMM2(0)
;day2=YYYYDDMM2(1)
;min2=YYYYDDMM2(2)
time2=SS2
;=================================================;
; On mets toutes les variables sur la grille 1 
;=================================================;


lsum1 = wrf_user_unstagger(lsum1_old,"X")
lsum1@long_name="LSUM"
lsum1@units="m/s"
ut1 = wrf_user_unstagger(ut1_old,"X")
ut1@long_name="Zonal wind"
ut1@units="m/s"

lsum2 = wrf_user_unstagger(lsum1_old,"X")
lsum2@long_name="LSUM"
lsum2@units="m/s"
ut2 = wrf_user_unstagger(ut1_old,"X")
ut2@long_name="Zonal wind"
ut2@units="m/s"

; Unstagger wt (from grid 4 to 1)
    wt1=new((/kmax,imax/),double)
    wt2=new((/kmax,imax/),double)

    do k=0,kmax-2
     wt1(k,:)=(wt1_old(k,:)+wt1_old(k+1,:))/2.
     wt2(k,:)=(wt2_old(k,:)+wt2_old(k+1,:))/2.
    end do
     wt1(kmax-1,:)=2*wt1_old(kmax-1,:)-wt1_old(kmax-2,:)
     wt2(kmax-1,:)=2*wt2_old(kmax-1,:)-wt2_old(kmax-2,:)

wt1@long_name="Vertical wind"
wt1@units="m/s"
wt2@long_name="Vertical wind"
wt2@units="m/s"

;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
    nzh=new(kmax,double)
    do k=0,kmax-2
     nzh(k)=(zhat(k)+zhat(k+1))/2.
    end do
     nzh(kmax-1)=2*zhat(kmax-1)-zhat(kmax-2)

; Create a (altitude des niveaux modèle)
    z=new(dimsizes(tht1),double)
    zcoef=new(imax,double)
    zcoef=1.-zs/nzh(kmax-1)

    do i=0,imax-1
       z(:,i) = nzh*zcoef(i)+zs(i)
    end do

   xconf=conform(tht1,xhat(0:imax-1),1)
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

  resmap@tmXBLabelStride = 2   ; to reduce the number of labels on xaxis
  resmap@lbLabelStride = 2.   ; to reduce the number of labels on labelbar


;pour prendre en compte l'orographie sur la coupe verticale.
  resmap@sfYArray        = z             ; 2D                  
  resmap@sfXArray        = xconf             ; 2D                  
  resmap@trGridType            = "TriangularMesh"
  

resmap@tiXAxisPosition="Left"
resmap@tiXAxisFontHeightF=0.015
;=================================================;
; TRACE
;=================================================;
  gsn_define_colormap(wks,"rainbow") ; Choose colormap

opts=resmap
opts@cnLevelSelectionMode = "ExplicitLevels"
cnLevels  = (/4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12/)
opts@cnLevels    = cnLevels
opts@tiXAxisString="time="+time1+"s"
plot_ut1 = gsn_csm_contour(wks,ut1,opts)
draw(plot_ut1)
frame(wks)
opts@tiXAxisString="time="+time2+"s"
plot_ut2 = gsn_csm_contour(wks,ut2,opts)
draw(plot_ut2)
frame(wks)
delete(opts)
delete(cnLevels)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 2.5
opts@cnMinLevelValF    = 285
opts@cnMaxLevelValF    = 330

opts@tiXAxisString="time="+time1+"s"

plot_tht1 = gsn_csm_contour(wks,tht1,opts)
draw(plot_tht1)
frame(wks)
opts@tiXAxisString="time="+time2+"s"
plot_tht2 = gsn_csm_contour(wks,tht2,opts)
draw(plot_tht2)
frame(wks)
delete(opts)


  resmap@lbLabelStride = 3.   ; to reduce the number of labels on labelbar
opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 0.04
opts@cnMinLevelValF    = -0.4
opts@cnMaxLevelValF    = 0.28
opts@tiXAxisString="time="+time1+"s"
plot_wt1 = gsn_csm_contour(wks,wt1,opts)
draw(plot_wt1)
frame(wks)
opts@tiXAxisString="time="+time2+"s"
plot_wt2 = gsn_csm_contour(wks,wt2,opts)
draw(plot_wt2)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 0.0002
opts@cnMinLevelValF    = 0.0002
opts@cnMaxLevelValF    = 0.0038
opts@tiXAxisString="time="+time1+"s"
;opts@lbLabelSride = 5.
plot_rvt1 = gsn_csm_contour(wks,rvt1,opts)
draw(plot_rvt1)
frame(wks)
opts@tiXAxisString="time="+time2+"s"
plot_rvt2 = gsn_csm_contour(wks,rvt2,opts)
draw(plot_rvt2)
frame(wks)
delete(opts)

  resmap@lbLabelStride = 2.   ; to reduce the number of labels on labelbar
opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 0.003
opts@cnMinLevelValF    = 0.028
opts@cnMaxLevelValF    = 0.078
opts@tiXAxisString="time="+time1+"s"
plot_cflu1 = gsn_csm_contour(wks,cflu1,opts)
draw(plot_cflu1)
frame(wks)

opts@tiXAxisString="time="+time2+"s"
plot_cflu2 = gsn_csm_contour(wks,cflu2,opts)
draw(plot_cflu2)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF    = 0.004
opts@cnMinLevelValF    = 0.004
opts@cnMaxLevelValF    = 0.072
opts@tiXAxisString="time="+time1+"s"
plot_cflw1 = gsn_csm_contour(wks,cflw1,opts)
draw(plot_cflw1)
frame(wks)
opts@tiXAxisString="time="+time2+"s"
plot_cflw2 = gsn_csm_contour(wks,cflw2,opts)
draw(plot_cflw2)
frame(wks)
delete(opts)


res=resmap
;===================================================
; profils verticaux
;===================================================
;====================================
 resmap@pmLegendDisplayMode    = "Always"             
 resmap@pmLegendWidthF         = 0.05                  
 resmap@pmLegendHeightF        = 0.15                  
 resmap@lgPerimOn              = False                 
 resmap@pmLegendSide           = "Top"                  
 resmap@pmLegendParallelPosF   = .1
 resmap@lgLabelFontHeightF     = .01                   
 resmap@pmLegendOrthogonalPosF = -0.3                
 resmap@tmXBLabelStride = 1   ; 

; delete ressources not needed for vertcal profiles
 delete(resmap@sfXArray)
 delete(resmap@sfYArray)
 delete(resmap@cnFillOn)
 delete(resmap@cnLinesOn)
 delete(resmap@gsnSpreadColors)
 delete(resmap@lbLabelStride)


res_w=resmap
res_t=resmap

res_w@gsnRightString="time="+time1+"s"
res_w@trXMinF=-0.3
res_w@trXMaxF=0.3
res_w@xyLineColors      = (/"red","green"/)          ; change line color
res_w@xyExplicitLegendLabels = (/"WT","THT-LSTHM"/)
res_w@tiXAxisString      = "WT"          ; title
res_w@tiXAxisFontHeightF = 0.015
res_t@pmLegendOrthogonalPosF = -0.4
res_t@pmLegendParallelPosF   = .15     
res_t@trXMinF=-0.8
res_t@trXMaxF=0.1
res_t@xyLineColors      = (/"green"/)          ; change line color
res_t@xyExplicitLegendLabels = (/"THT-LSTHM"/)
res_t@gsnLeftString      = "THT-LSTHM"          ; title
res_t@gsnLeftStringFontHeightF = 0.015

plot  = gsn_csm_x2y (wks,wt1(:,127),tht1(:,127)-lsthm1(:,127),z(:,127),res_w,res_t)            
draw(plot)
frame(wks)

res_t@trXMinF=-1
res_t@trXMaxF=0.2
res_w@gsnRightString="time="+time2+"s"
plot  = gsn_csm_x2y (wks,wt2(:,127),tht2(:,127)-lsthm2(:,127),z(:,127),res_w,res_t)            
draw(plot)
frame(wks)

res_u=resmap
res_lsu=resmap

res_w@xyExplicitLegendLabels = (/"WT","UT-LSUM"/)
res_lsu@pmLegendOrthogonalPosF = -0.4
res_lsu@pmLegendParallelPosF   = .15     
res_lsu@trXMinF=-0.2
res_lsu@trXMaxF=0.8
res_lsu@xyLineColors      = (/"green"/)          ; change line color
res_lsu@xyExplicitLegendLabels = (/"UT-LSUM"/)
res_lsu@gsnLeftString      = "UT-LSUM"          ; title
res_lsu@tiMainFontHeightF = 0.015
res_w@gsnRightString="time="+time1+"s"

plot  = gsn_csm_x2y (wks,wt1(:,127),ut1(:,127)-lsum1(:,127),z(:,127),res_w,res_lsu)            
draw(plot)
frame(wks)
res_lsu@trXMinF=-2
res_lsu@trXMaxF=2
res_w@gsnRightString="time="+time2+"s"

plot  = gsn_csm_x2y (wks,wt2(:,127),ut2(:,127)-lsum2(:,127),z(:,127),res_w,res_lsu)            
draw(plot)
frame(wks)


end
