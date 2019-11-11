;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
; ================================================;
begin
;=================================================;
; open file and read in data
;=================================================;
  a = addfile("ICART.1.SEG01.001dg.nc", "r")
  a2 = addfile("ICART.1.SEG01.002dg.nc", "r")

;=================================================;
; Get informations on variable sizes
; dims are dims-2 to remove non-physical values
;=================================================;
  jphext  = a->JPHEXT
  mdims = getfilevardimsizes(a,"THT") ; get dimension sizes
  nd = dimsizes(mdims)
  imax=mdims(nd-1)-2*jphext
  jmax=mdims(nd-2)-2*jphext
  kmax=mdims(nd-3)-2

;-------------------------------------------------;
; Read data.
;-------------------------------------------------;

; Liste de toutes les réactions impliquant O3
  o3_list_t1=a->O3_CHREACLIST(0,:)

; Tableau 4D (reac,Z,Y,X) regroupant les termes de prod. et destr. de O3
; Niveau k=13 approx. 1250m
  o3_budget_t1=a->O3_BUDGET(0,:,13,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  printVarSummary(o3_budget_t1)
  o3_bud_t1=dim_avg_n(o3_budget_t1,(/1,2/))
   printVarSummary(o3_bud_t1)
 
  o3_bud_t1=o3_bud_t1*1e9*3600
  o3_bud_t1@units="ppbv/h"

; Liste de toutes les réactions impliquant CO
  co_list_t1=a->CO_CHREACLIST(0,:)

; Tableau 4D (reac,Z,Y,X) regroupant les termes de prod. et destr. de CO
; Niveau k=13 approx. 1250m
  co_budget_t1=a->CO_BUDGET(0,:,13,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  co_bud_t1=dim_avg_n(co_budget_t1,(/1,2/))
  co_bud_t1=co_bud_t1*1e9*3600
  co_bud_t1@units="ppbv/h"

; Liste de toutes les réactions impliquant O3
  o3_list_t2=a2->O3_CHREACLIST(0,:)

; Tableau 4D (reac,Z,Y,X) regroupant les termes de prod. et destr. de O3
; Niveau k=13 approx. 1250m
  o3_budget_t2=a2->O3_BUDGET(0,:,13,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  o3_bud_t2=dim_avg_n(o3_budget_t2,(/1,2/))
  o3_bud_t2=o3_bud_t2*1e9*3600
  o3_bud_t2@units="ppbv/h"

; Liste de toutes les réactions impliquant CO
  co_list_t2=a2->CO_CHREACLIST(0,:)

; Tableau 4D (reac,Z,Y,X) regroupant les termes de prod. et destr. de CO
; Niveau k=13 approx. 1250m
  co_budget_t2=a2->CO_BUDGET(0,:,13,jphext:jmax+jphext-1,jphext:imax+jphext-1)
  co_bud_t2=dim_avg_n(co_budget_t2,(/1,2/))
  co_bud_t2=co_bud_t2*1e9*3600
  co_bud_t2@units="ppbv/h"

;=================================================;
; PLOT
;=================================================;
; interpolation des champs a 1250 m

  figname ="zsection_1250_bud"
  wks  = gsn_open_wks("png",figname)   ; open a ncgm file
  gsn_define_colormap(wks,"WhBlGrYeRe") ; Choose colormap

  res                 = True          
  res@gsnDraw   		= False         ; don't draw yet
  res@gsnFrame  		= False         ; don't advance frame yet

; X-axis title (tiY)                              
  res@tiXAxisFontHeightF = 0.018  ; font height
  res@tiXAxisFont        = 21     ; font index
  res@tiXAxisString      = "Chemical reactions"  ; string to use as the X-Axis title

; Y-axis title (tiY)
  res@tiYAxisFontHeightF = 0.018  ; font height
  res@tiYAxisFont        = 21     ; font index

; Bar plot
  res@gsnXYBarChart         = True                 ; turn on bar chat
  res@gsnYRefLine           = 0.              ; reference line
  res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
  res@xyCurveDrawOrder      = "PreDraw"  
  res@tmYLFormat            = "0@*+^sg"
  res@tmYLPrecision         = 2
  res@tmXBOn=False

  txres                     = True                 ; text mods desired
  txres@txFontHeightF       = 0.018                 
  txres@txAngleF            = 90                  ; text angle
  txres@txJust              = "TopCenter"         ; puts text on top of bars

  res@gsnCenterString="heure=19"

; plot ozone production
  res@tiYAxisString      = "ozone budget (ppbv/h)"  ; string to use as the Y-Axis title
  x=ispan(0,dimsizes(o3_bud_t1)-1,1)
  plot_o3_t1 = gsn_csm_xy(wks,x,o3_bud_t1(:),res)
  text=gsn_add_text(wks,plot_o3_t1,tostring(o3_list_t1(:)),x,o3_bud_t1(:),txres) ; add label
  draw(plot_o3_t1)
  frame(wks)

; plot ozone production
  res@tiYAxisString      = "carbon monoxide budget (ppbv/h)"  ; string to use as the Y-Axis title
  x:=ispan(0,dimsizes(co_bud_t1)-1,1)
  plot_co_t1 = gsn_csm_xy(wks,x,co_bud_t1(:),res)
  text:=gsn_add_text(wks,plot_co_t1,tostring(co_list_t1(:)),x,co_bud_t1(:),txres) ; add label
  draw(plot_co_t1)
  frame(wks)

  res@gsnCenterString="heure=20"

; plot ozone production
  res@tiYAxisString      = "ozone budget (ppbv/h)"  ; string to use as the Y-Axis title
  x:=ispan(0,dimsizes(o3_bud_t2)-1,1)
  plot_o3_t2 = gsn_csm_xy(wks,x,o3_bud_t2(:),res)
  text:=gsn_add_text(wks,plot_o3_t1,tostring(o3_list_t2(:)),x,o3_bud_t2(:),txres) ; add label
  draw(plot_o3_t2)
  frame(wks)

; plot ozone production
  res@tiYAxisString      = "carbon monoxide budget (ppbv/h)"  ; string to use as the Y-Axis title
  x:=ispan(0,dimsizes(co_bud_t2)-1,1)
  plot_co_t2 = gsn_csm_xy(wks,x,co_bud_t2(:),res)
  text:=gsn_add_text(wks,plot_co_t2,tostring(co_list_t2(:)),x,co_bud_t2(:),txres) ; add label
  draw(plot_co_t2)
  frame(wks)



;;;;;;;;;;;;;;;;;;;;;;;;

end



