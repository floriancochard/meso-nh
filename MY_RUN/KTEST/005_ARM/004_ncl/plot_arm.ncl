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
  fichier1 = addfile("ARM__.1.CEN4T.000.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_ARM")
;=================================================;
; lecture des différents champs
;=================================================;
kmax=100
zhat  = fichier1->ZHAT(:) ; ZHAT
; Unstagger zhat (from grid 4 to 1)
    nzh=new(kmax,double)
    do k=0,kmax-2
     nzh(k)=(zhat(k)+zhat(k+1))/2.
    end do
     nzh(kmax-1)=2*zhat(kmax-1)-zhat(kmax-2)
TIME=new((/180/),double)
TIME(0)=300
do j=1,179
 TIME(j)=TIME(j-1)+300
end do



vname="SBG_TKE___PROC1"
SBGTKE = fichier1->$vname$(0,:,:,0,0) 
SBGTKE@long_name="SBG_TKE"
SBGTKE!1="nzh"
SBGTKE!0="TIME"
SBG_TKE=transpose(SBGTKE)

vname="MEAN_TH___PROC1"
MEANTH = fichier1->$vname$(0,:,:,0,0) ; MEAN_TH
MEANTH@long_name="MEAN_TH"
MEANTH@units="K"
MEANTH!1="nzh"
MEANTH!0="TIME"
MEAN_TH=transpose(MEANTH)


vname="MEAN_U___PROC1"
MEANU = fichier1->$vname$(0,:,:,0,0) 
MEANU!1="nzh"
MEANU!0="TIME"
MEAN_U=transpose(MEANU)
MEAN_U@long_name="MEAN_U"

vname="MEAN_V___PROC1"
MEANV = fichier1->$vname$(0,:,:,0,0) 
MEANV!1="nzh"
MEANV!0="TIME"
MEAN_V=transpose(MEANV)
MEAN_V@long_name="MEAN_V"

vname="MEAN_RC___PROC1"
MEANRC = fichier1->$vname$(0,:,:,0,0) 
MEANRC!1="nzh"
MEANRC!0="TIME"
MEAN_RC=transpose(MEANRC)
MEAN_RC@long_name="MEAN_RC"

vname="MEAN_RR___PROC1"
MEANRR = fichier1->$vname$(0,:,:,0,0) 
MEANRR!1="nzh"
MEANRR!0="TIME"
MEAN_RR=transpose(MEANRR)
MEAN_RR@long_name="MEAN_RR"

vname="SBG_WTHL___PROC1"
SBGWTHL = fichier1->$vname$(0,:,:,0,0) 
SBGWTHL!1="nzh"
SBGWTHL!0="TIME"
SBG_WTHL=transpose(SBGWTHL)
SBG_WTHL@long_name="SBG_WTHL"

vname="SBG_WRT___PROC1"
SBGWRT = fichier1->$vname$(0,:,:,0,0) 
SBGWRT!1="nzh"
SBGWRT!0="TIME"
SBG_WRT=transpose(SBGWRT)
SBG_WRT@long_name="SBG_WRT"


vname="THLUP_MF___PROC1"
 THLUPMF= fichier1->$vname$(0,:,:,0,0) 
THLUPMF!1="nzh"
THLUPMF!0="TIME"
THLUP_MF=transpose(THLUPMF)
THLUP_MF@long_name="THLUP_MF"

vname="RTUP_MF___PROC1"
RTUPMF = fichier1->$vname$(0,:,:,0,0) 
RTUPMF!1="nzh"
RTUPMF!0="TIME"
RTUP_MF=transpose(RTUPMF)
RTUP_MF@long_name="RTUP_MF"

vname="RVUP_MF___PROC1"
RVUPMF = fichier1->$vname$(0,:,:,0,0) 
RVUPMF!1="nzh"
RVUPMF!0="TIME"
RVUP_MF=transpose(RVUPMF)
RVUP_MF@long_name="RVUP_MF"

vname="RCUP_MF___PROC1"
RCUPMF = fichier1->$vname$(0,:,:,0,0) 
RCUPMF!1="nzh"
RCUPMF!0="TIME"
RCUP_MF=transpose(RCUPMF)
RCUP_MF@long_name="RCUP_MF"

vname="RIUP_MF___PROC1"
RIUPMF = fichier1->$vname$(0,:,:,0,0) 
RIUPMF!1="nzh"
RIUPMF!0="TIME"
RIUP_MF=transpose(RIUPMF)
RIUP_MF@long_name="RIUP_MF"

vname="WUP_MF___PROC1"
WUPMF = fichier1->$vname$(0,:,:,0,0) 
WUPMF!1="nzh"
WUPMF!0="TIME"
WUP_MF=transpose(WUPMF)
WUP_MF@long_name="WUP_MF"

vname="MAFLX_MF___PROC1"
MAFLXMF = fichier1->$vname$(0,:,:,0,0) 
MAFLXMF!1="nzh"
MAFLXMF!0="TIME"
MAFLX_MF=transpose(MAFLXMF)
MAFLX_MF@long_name="MAFLX_MF"

vname="DETR_MF___PROC1"
DETRMF = fichier1->$vname$(0,:,:,0,0) 
DETRMF!1="nzh"
DETRMF!0="TIME"
DETR_MF=transpose(DETRMF)
DETR_MF@long_name="DETR_MF"

vname="ENTR_MF___PROC1"
ENTRMF = fichier1->$vname$(0,:,:,0,0) 
ENTRMF!1="nzh"
ENTRMF!0="TIME"
ENTR_MF=transpose(ENTRMF)
ENTR_MF@long_name="ENTR_MF"

vname="FRCUP_MF___PROC1"
FRCUPMF = fichier1->$vname$(0,:,:,0,0) 
FRCUPMF!1="nzh"
FRCUPMF!0="TIME"
FRCUP_MF=transpose(FRCUPMF)
FRCUP_MF@long_name="FRCUP_MF"

vname="THVUP_MF___PROC1"
THVUPMF = fichier1->$vname$(0,:,:,0,0) 
THVUPMF!1="nzh"
THVUPMF!0="TIME"
THVUP_MF=transpose(THVUPMF)
THVUP_MF@long_name="THVUP_MF"

vname="WTHL_MF___PROC1"
WTHLMF = fichier1->$vname$(0,:,:,0,0) 
WTHLMF!1="nzh"
WTHLMF!0="TIME"
WTHL_MF=transpose(WTHLMF)
WTHL_MF@long_name="WTHL_MF"

vname="WRT_MF___PROC1"
WRTMF = fichier1->$vname$(0,:,:,0,0) 
WRTMF!1="nzh"
WRTMF!0="TIME"
WRT_MF=transpose(WRTMF)
WRT_MF@long_name="WRT_MF"


vname="WTHV_MF___PROC1"
WTHVMF = fichier1->$vname$(0,:,:,0,0) 
WTHVMF!1="nzh"
WTHVMF!0="TIME"
WTHV_MF=transpose(WTHVMF)
WTHV_MF@long_name="WTHV_MF"


vname="WU_MF___PROC1"
WUMF = fichier1->$vname$(0,:,:,0,0) 
WUMF!1="nzh"
WUMF!0="TIME"
WU_MF=transpose(WUMF)
WU_MF@long_name="WU_MF"


vname="WV_MF___PROC1"
WVMF = fichier1->$vname$(0,:,:,0,0) 
WVMF!1="nzh"
WVMF!0="TIME"
WV_MF=transpose(WVMF)
WV_MF@long_name="WV_MF"

;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
    nzh=new(kmax,double)
    do k=0,kmax-2
     nzh(k)=(zhat(k)+zhat(k+1))/2.
    end do
     nzh(kmax-1)=2*zhat(kmax-1)-zhat(kmax-2)
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
  ;resmap@tmXBLabelStride = 2   ; to reduce the number of labels on xaxis
  ;resmap@lbLabelStride = 2.   ; to reduce the number of labels on labelbar
  resmap@sfYArray        = nzh(0:76)                            
  resmap@sfXArray        = TIME               
  resmap@tmXBTickSpacingF  = 10800.
  resmap@tmXBMode          = "Explicit"
  resmap@tmXBValues = (/10800,21600,32400,43200,54000/)
  resmap@tmXBLabels =(/3,6,9,12,15/)
resmap@tiXAxisPosition="Left"
resmap@tiXAxisFontHeightF=0.015
;=================================================;
; TRACE
;=================================================;
  gsn_define_colormap(wks,"rainbow") ; Choose colormap
       colors = gsn_retrieve_colormap(wks)     ; retrieve color map for editing
       colors(2,:) = (/ 1, 1, 1 /)     ; replace the first color with white color (les deux premières sont background /foreground donc c'est bien la 2 qu'il faut changer)
       gsn_define_colormap(wks,colors)         ; redefine colormap to workstation, color map now includes a gray



opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.25
opts@cnMinLevelValF=0
opts@cnMaxLevelValF=2.5

plot = gsn_csm_contour(wks,SBG_TKE(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=1
opts@cnMinLevelValF=298
opts@cnMaxLevelValF=316

plot = gsn_csm_contour(wks,MEAN_TH(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)


opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.5
opts@cnMinLevelValF=3
opts@cnMaxLevelValF=11

plot = gsn_csm_contour(wks,MEAN_U(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.2
opts@cnMinLevelValF=0.
opts@cnMaxLevelValF=3.2

plot = gsn_csm_contour(wks,MEAN_V(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.00001
opts@cnMinLevelValF=0
opts@cnMaxLevelValF=0.0002

plot = gsn_csm_contour(wks,MEAN_RC(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)


opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.000001
opts@cnMinLevelValF=0
opts@cnMaxLevelValF=0.00002

plot = gsn_csm_contour(wks,MEAN_RR(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)


opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.02
opts@cnMinLevelValF=-0.3
opts@cnMaxLevelValF=0.12

plot = gsn_csm_contour(wks,SBG_WTHL(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.00005
opts@cnMinLevelValF=-0.0001
opts@cnMaxLevelValF=0.0004

plot = gsn_csm_contour(wks,SBG_WRT(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)


opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=1
opts@cnMinLevelValF=301
opts@cnMaxLevelValF=318

plot = gsn_csm_contour(wks,THLUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.001
opts@cnMinLevelValF=0.0035
opts@cnMaxLevelValF=0.017
plot = gsn_csm_contour(wks,RTUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.001
opts@cnMinLevelValF=0.001
opts@cnMaxLevelValF=0.017

plot = gsn_csm_contour(wks,RVUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.00025
opts@cnMinLevelValF=0.00025
opts@cnMaxLevelValF=0.00425

plot = gsn_csm_contour(wks,RCUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "AutomaticLevels"
plot = gsn_csm_contour(wks,RIUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.25
opts@cnMinLevelValF=0.250
opts@cnMaxLevelValF=5.5

plot = gsn_csm_contour(wks,WUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.02
opts@cnMinLevelValF=0.02
opts@cnMaxLevelValF=0.32

plot = gsn_csm_contour(wks,MAFLX_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.025
opts@cnMinLevelValF=0.025
opts@cnMaxLevelValF=0.45

plot = gsn_csm_contour(wks,DETR_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.00025
opts@cnMinLevelValF=0.00025
opts@cnMaxLevelValF=0.005

plot = gsn_csm_contour(wks,ENTR_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.01
opts@cnMinLevelValF=0.01
opts@cnMaxLevelValF=0.16

plot = gsn_csm_contour(wks,FRCUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=1
opts@cnMinLevelValF=303
opts@cnMaxLevelValF=319

plot = gsn_csm_contour(wks,THVUP_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.005
opts@cnMinLevelValF=-0.075
opts@cnMaxLevelValF=0.055

plot = gsn_csm_contour(wks,WTHL_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.00001
opts@cnMinLevelValF=0.00001
opts@cnMaxLevelValF=0.00019

plot = gsn_csm_contour(wks,WRT_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.005
opts@cnMinLevelValF=-0.016
opts@cnMaxLevelValF=0.075
plot = gsn_csm_contour(wks,WTHV_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.005
opts@cnMinLevelValF=-0.120
opts@cnMaxLevelValF=-0.005

plot = gsn_csm_contour(wks,WU_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

opts=resmap
opts@cnLevelSelectionMode = "ManualLevels"
opts@cnLevelSpacingF=0.004
opts@cnMinLevelValF=-0.024
opts@cnMaxLevelValF=0.048
plot = gsn_csm_contour(wks,WV_MF(0:76,:),opts)
draw(plot)
frame(wks)
delete(opts)

end
