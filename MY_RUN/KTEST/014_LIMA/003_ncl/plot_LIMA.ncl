;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; ================================================;
;================================================;
begin
;=================================================;
; open file and read in data
;=================================================;
  fichier1 = addfile("XPREF.1.SEG01.000.nc", "r")
;==================================================;
; Open the workstation
;==================================================;
  type = "png"
  wks = gsn_open_wks(type,"visu_LIMA")
;=================================================;
; lecture des différents champs
;=================================================;
vname = "RI___0001___PROC3"
RI___0001_P3= fichier1->$vname$(0,0,:,0,:)
RI___0001_P3@long_name="ice water content"
RI___0001_P3@units="g/kg"

vname2 = "SV0050001___PROC3"
NI___0001_P3= fichier1->$vname2$(0,0,:,0,:)
NI___0001_P3@long_name="ice concentration"
NI___0001_P3@units="/kg"

vname3 = "RS___0001___PROC3"
RS___0001_P3= fichier1->$vname3$(0,0,:,0,:)
RS___0001_P3@long_name="snow water content"
RS___0001_P3@units="g/kg"

vname4 = "RG___0001___PROC3"
RG___0001_P3= fichier1->$vname4$(0,0,:,0,:)
RG___0001_P3@long_name="graupel water content"
RG___0001_P3@units="g/kg"

vname5 = "SV0060001___PROC3"
N_IFN_FREE_P3= fichier1->$vname5$(0,0,:,0,:)
N_IFN_FREE_P3@long_name="concentration of free IFN"
N_IFN_FREE_P3@units="/kg"

vname6 = "SV0070001___PROC3"
N_IFN_NUCL_P3= fichier1->$vname6$(0,0,:,0,:)
N_IFN_NUCL_P3@long_name="concentration of nucleated IFN"
N_IFN_NUCL_P3@units="/kg"


zhat  = fichier1->ZHAT(1:50) 
xhat  = fichier1->XHAT(1:180) 
zs  = fichier1->ZS(0,1:180) ; ZS

xconf=conform(RI___0001_P3,xhat(0:179),1)

;=================================================;
; Altitude des niveaux modèles
;=================================================;
; Unstagger zhat (from grid 4 to 1)
    nzh=new(50,double)
    do k=0,48
     nzh(k)=(zhat(k)+zhat(k+1))/2.
    end do
     nzh(49)=2*zhat(49)-zhat(48)

; Create z (altitude des niveaux modèle)
    z=new((/50,180/),double)
    zcoef=new(180,double)
    zcoef=1.-zs/nzh(49)

    do i=0,179
       z(:,i) = nzh*zcoef(i)+zs(i)
    end do
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

  resmap@tiXAxisPosition="Left"
  resmap@tiXAxisFontHeightF=0.015
;  resmap@gsnSpreadColorStart  = 0     ; force la première couleur en blanc(= couleur 0 de la palette)
  resmap@sfYArray        = z             ; 2D                  
  resmap@sfXArray        = xconf             ; 2D  
  resmap@trGridType            = "TriangularMesh"  
;=================================================;
; TRACE
;=================================================;
  gsn_define_colormap(wks,"rainbow") ; Choose colormap
  res=resmap       
  res@cnLevelSelectionMode = "AutomaticLevels"
  plot = gsn_csm_contour(wks,RI___0001_P3(:,:),res)
  draw(plot)
  frame(wks)

  plot2 = gsn_csm_contour(wks,NI___0001_P3(:,:),res)
  draw(plot2)
  frame(wks)

  plot3 = gsn_csm_contour(wks,RS___0001_P3(:,:),res)
  draw(plot3)
  frame(wks)

  plot4 = gsn_csm_contour(wks,RG___0001_P3(:,:),res)
  draw(plot4)
  frame(wks)
  plot5 = gsn_csm_contour(wks,N_IFN_FREE_P3(:,:),res)
  draw(plot5)
  frame(wks)
  plot6 = gsn_csm_contour(wks,N_IFN_NUCL_P3(:,:),res)
  draw(plot6)
  frame(wks)
 
end
