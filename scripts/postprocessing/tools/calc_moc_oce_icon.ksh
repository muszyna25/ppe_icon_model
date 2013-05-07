#!/bin/ksh

#  Calculate meridional overturning stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 07/2012 - 04/2013
#
#  Method:
#   - calculation of global/Atlantic/Pacific+Indian MOC is done in the model
#   - written in extra format
#   - simple splitting of MOC result file and plotting using cdo and ncl
#   - options: give masked file, which can be computed from the input like this:
#     cdo -setmisstoc,0.0 -div timmean_moc_0185-0190.ext timmean_moc_0185-0190.ext mocMask.nc

ifile=moc.bliz.r9368.2420-2470ym
ifile=moc.loc.r9558.2046.2y
ifile=moc.loc.r9513.2042.2y
ifile=moc.loc.r9558.2061.2y
ifile=$1
maskfile=$2

# check intput:
echo "Input file is '$ifile'"

basename=$(basename $ifile .ext) # ext files expected

if [ -z "$maskfile" ]; then
  cdo -f nc selvar,var777 $ifile scr_moc_glb.nc 
  cdo -f nc selvar,var778 $ifile scr_moc_atl.nc
  cdo -f nc selvar,var779 $ifile scr_moc_pac.nc
else
  cdo -f nc selvar,var777 $ifile scr_moc_glb.nc 
  cdo -div scr_moc_glb.nc -selname,var777 $maskfile tmp.nc
  mv tmp.nc scr_moc_glb.nc

  cdo -f nc selvar,var778 $ifile scr_moc_atl.nc
  cdo -div scr_moc_atl.nc -selname,var778 $maskfile tmp.nc
  mv tmp.nc scr_moc_atl.nc

  cdo -f nc selvar,var779 $ifile scr_moc_pac.nc
  cdo -div scr_moc_pac.nc -selname,var779 $maskfile tmp.nc
  mv tmp.nc scr_moc_pac.nc
fi

# run ncl-script

cat >scr_plot_moc_my.ncl << EOF
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************

begin

; if (isvar("infilename")) then
;   prefix=infilename + "_"
; else
; prefix=""
; end if
         
 f = addfile("scr_moc_atl.nc","r")
 g = addfile("scr_moc_glb.nc","r")
 h = addfile("scr_moc_pac.nc","r")

 moc_atl     = f->var778(0,:,60:179)
 moc_glo     = g->var777(0,:,0:179)
 moc_pac_ind = h->var779(0,:,40:159)

 moc_atl     = moc_atl/1000000000
 moc_glo     = moc_glo/1000000000
 moc_pac_ind = moc_pac_ind/1000000000

;************************************************************************
;  more TBD:
;   - 60S instead of -60 -> needs correct grid description
;************************************************************************

  wks = gsn_open_wks("ps","$basename")
; wks = gsn_open_wks("pdf","$basename")

  gsn_define_colormap(wks,"BlWhRe")            ; colormap

  res                      = True
; res@cnFillMode           = "RasterFill"
  res@tmXBMode 		   = "Explicit"
; res@gsnPaperOrientation  = "Portrait"
  res@gsnPaperOrientation  = "Landscape"
  res@gsnMaximize          = True              ; important

  res@gsnYAxisIrregular2Linear = True          ; True: axis-scaling is linear in Y, False: following resolution

  res@vpWidthF             = 0.8               ; change aspect ratio of plot
  res@vpHeightF            = 0.45              ; change aspect ratio of plot
  res@cnFillOn             = True              ; color for contours
; res@cnLinesOn            = False             ; contour lines
  res@lbLabelBarOn         = False             ; color label bar
  res@cnInfoLabelOn        = False             ; contour info label
  res@gsnSpreadColors      = True              ; subset of the whole colormap
; res@gsnSpreadColorStart  =  25               ; start color
; res@gsnSpreadColorEnd    = -25               ; end color
; res@cnFillMode           = "CellFill"        ; filled cells


; choice: range of levels
; res@cnLevelSelectionMode  = "ManualLevels"    ; set manual contour levels
; res@cnMinLevelValF        = -30               ; set min contour level
; res@cnMaxLevelValF        = 30                ; set max contour level
; res@cnLevelSpacingF       = 2                 ; set contour spacing

; choice: explicit levels
  res@cnLevelSelectionMode  = "ExplicitLevels"  ; set explicit contour levels
; res@cnLevels              = (/ -100, -20, -10, -5,0,5,10,20,100 /)
; res@cnLevels              = (/ -200, -150, -100, -50, -30, -20, -15, -10, -7, -5, -2, 0, 2, 5, 7, 10, 15, 20, 50, 100, 150, 200 /)
; res@cnLevels              = (/ -100, -50, -30, -25, -20, -16, -13, -10, -7, -5, -3, -1, 0, 1, 3, 5, 7, 10, 13, 16, 20, 25, 30, 50, 100 /)
; res@cnLevels              = (/ -100, -80, -60, -50, -40, -35, -30, -25, -22.5, -20, -18, -16, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, \
;                                1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22.5, 25, 30, 35, 40, 50, 60, 80, 100 /)
; res@cnLevels              = (/ -50, -40, -35, -30, -25, -22.5, -20, -18, -16, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, \
;                                1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22.5, 25, 30, 35, 40, 50 /)
  res@cnLevels              = (/ -27.5, -25, -22.5, -20, -18, -16, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, \
                                 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22.5, 25, 27.5 /)

  res@gsnContourZeroLineThicknessF   = 3.0      ; factor for thickness of zero contour
  res@gsnContourLineThicknessesScale = 2.0      ; factor for other lines

  res@cnLineLabelsOn        = True              ; contour line labels
; res@cnLineLabelsOn        = False             ; contour line labels
  res@cnLineLabelFontHeightF= 0.013             ; contour line label font
  res@cnLineLabelInterval   = 2                 ; line label interval (default=2)
  res@cnLabelMasking        = True		; masking contour lines below line labels
  res@cnLineLabelBackgroundColor = "transparent"
  res@cnLineLabelBackgroundColor = "white"

  res@lbLabelFontHeightF    = 0.014             ; label bar font
  res@lbLabelStride         = 2                 ; labelling interval at color bar
  res@lbLabelPosition       = "Left"  
; res@lbOrientation         = "horizontal"
; res@lbOrientation	    = "vertical" 
; res@pmLabelBarOrthogonalPosF = 0.1            ; move label bar closer

 
  res@cnMissingValFillColor        = "gray30"    ; not defined?
  res@cnMissingValPerimOn          = True
  res@cnMissingValFillPattern      = -1          ; set the missing value fill pattern
  res@cnMissingValFillScaleF       = 0.9         ; increase the density of the fill pattern (default   = 1.0)
  res@cnMissingValPerimColor       = "black"     ; change the missing value perimeter to black
  res@cnMissingValPerimDashPattern = 1           ; set the dash pattern of the missing value perimeter to 1
  res@cnMissingValPerimThicknessF  = 3.0         ; factor for thickness of the missing value perimeter 3X

  res@gsnRightString	       = "[Sv]"
  res@gsnLeftString	       = "$basename"     ; filename string
  res@trYReverse               = True            ; reverse the Y-axis
  res@tiXAxisString            = "latitude"
  res@tiYAxisString            = "depth [m]"

; *****************************************
; Plotting experiments
; *****************************************

  res@tmXBValues 	       = ispan(0,180,30)
  res@tmXBLabels 	       = ispan(-90,90,30)
; res@trXMinF                  = -80           ; set minimum X-axis value
; res@trXMaxF                  =  85           ; set maximum X-axis value
; res@tiMainString	       = "MOC Global "  ; + infilename
; res@tiMainString	       = "MOC Global - $basename "  ; + infilename
  res@tiMainString	       = "MOC Global"

  plot = gsn_csm_contour(wks,moc_glo,res)

; wks = gsn_open_wks("ps","$basename.atl")   ;  new resources are not complete for new filename
  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,120,30)
  res@tmXBLabels 	       = ispan(-30,90,30)
  res@tiMainString	       = "MOC Atlantic"
  plot = gsn_csm_contour(wks,moc_atl,res)

  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,110,30)
  res@tmXBLabels 	       = ispan(-50,70,30)
  res@tiMainString	       = "MOC Pac+Ind"
  plot = gsn_csm_contour(wks,moc_pac_ind,res)
end
EOF


ncl scr_plot_moc_my.ncl
rm scr_plot_moc_my.ncl
rm scr_moc_???.nc





