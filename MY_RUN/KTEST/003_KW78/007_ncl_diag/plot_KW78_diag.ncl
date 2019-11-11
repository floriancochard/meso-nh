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
  fichier1 = addfile("KWRAI.1.SEG01.004.nc", "r")
  fichier2 = addfile("KWRAI.1.SEG01.004dia.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_KW78_diag")
;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = fichier2->JPHEXT
  mdims = getfilevardimsizes(fichier2,"CLUSTERID") 
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
  
  idt1= fichier2->CLUSTERID(0,1:kmax,jphext:jmax+jphext-1,jphext:imax+jphext-1) ; CLUSTERID
  idt1@long_name="Cluster Identity Number"
  idt1@units=""

  YYYYDDMM1= fichier1->DTCUR
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
  resmap@gsnSpreadColors= True       
  resmap@tiYAxisString =" "
  resmap@cnFillOn= True 
  resmap@cnLinesOn= False 
  resmap@lbLabelStride = 12.   ; to reduce the number of labels on labelbar
;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
  nzh=new(kmax,double)
  do k=0,kmax-1
    nzh(k)=(zhat(k)+zhat(k+1))/2.
  end do

; Create alt (altitude des niveaux modèle)
  alt=new(dimsizes(ut1),double)
  zcoef=new((/jmax,imax/),double)
  zcoef=1.-zs/nzh(kmax-1)
  do i=0,imax-1
    do j=0,jmax-1
      alt(:,j,i) = nzh*zcoef(j,i)+zs(j,i)
    end do
  end do

  
;=======================================
; coupe verticale clustering
;=======================================

  gsn_define_colormap(wks,"BlRe")

  opt_cv=resmap
  opt_cv@gsnScalarContour=True
  opt_cv@vcMinDistanceF = 0.03 ; espacement des fleches
  opt_cv@vcRefLengthF=0.1      ; taille des flèches
  opt_cv@cnLevelSpacingF    = 0.0005
  opt_cv@cnMinLevelValF    = 0
  opt_cv@cnLevelSelectionMode = "ExplicitLevels"
  opt_cv@cnLevels         = (/1/)
  opt_cv@cnFillColors     = (/0,84/)
  opt_cv@tiMainString = "Identified Updraft as cluster where w > 10 m/s"
  
  plot_idt=gsn_csm_vector_scalar(wks,ut1(:,13,:),wt1(:,13,:),idt1(:,13,:),opt_cv)

 draw(plot_idt)
 frame(wks)
  

end
