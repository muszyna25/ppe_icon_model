#!/bin/ksh

#  Calculate meridional overturning stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 07/2012
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
                                                                                                                                                                                                          
;if (isvar("infilename")) then
;  prefix=infilename + "_"
;else
;prefix=""
;end if                                                                                                                                                                                                         
         
 f = addfile("scr_moc_atl.nc","r")
 g = addfile("scr_moc_glb.nc","r")
 h = addfile("scr_moc_pac.nc","r")

 moc_atl = f->var778(0,:,60:179)
 moc_glo  = g->var777(0,:,0:179)
 moc_pac_ind = h->var779(0,:,40:159)

 moc_atl = moc_atl/1000000000
 moc_glo = moc_glo/1000000000
 moc_pac_ind = moc_pac_ind/1000000000

;************************************************************************
;************************************************************************
  wks = gsn_open_wks("ps","$basename")
; wks = gsn_open_wks("pdf","$basename")

  gsn_define_colormap(wks,"BlWhRe")

  res                      = True
; res@cnFillMode             = "RasterFill"
  res@tmXBMode 		   = "Explicit"
; res@gsnPaperOrientation  = "Portrait"
  res@gsnPaperOrientation  = "Landscape"
  res@gsnMaximize          = True 

  res@vpWidthF             = 0.8           ; change aspect ratio of plot
  res@vpHeightF            = 0.5
  res@cnFillOn             = True              ; color for contours
; res@cnLinesOn            = False             ; contour lines
  res@gsnSpreadColors      = True              ; full color map
; res@cnFillMode           = "CellFill"

  res@cnLevelSelectionMode        = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF              = -19               ; set min contour level
  res@cnMaxLevelValF              = 19                ; set max contour level
  res@cnLevelSpacingF             = 2                 ; set contour spacing
; res@cnLevelSelectionMode        = "ExplicitLevels"  ; set explicit contour levels
; res@cnLevels                    = (/ -250, -150, -75, -50, -30, -20, -15, -10, -2, 0, 2, 10, 15, 20, 30, 50, 75, 150, 250 /)
  res@cnLevels                    = (/ -150, -100, -50, -30, -20, -15, -10, -7, -5, -2, 0, 2, 5, 7, 10, 15, 20, 30, 50, 100, 150 /)
; res@cnLevels                    = (/ -200, -150, -100, -50, -30, -20, -15, -10, -7, -5, -2, 0, 2, 5, 7, 10, 15, 20, 50, 100, 150, 200 /)
; res@cnLineLabelsOn              = False             ; contour line labels
; res@cnLevelSelectionMode        = "ExplicitLevels"  ; set explicit contour levels
; res@cnLevels                    = (/ -100, -20, -10, -5,0,5,10,20,100 /)
  res@cnLineLabelsOn              = True              ; contour line labels
  res@lbLabelStride        	  = 1
  res@lbLabelFontHeightF	  = 0.015
; res@lbOrientation               = "Vertical"
 
  res@cnLabelMasking             = True		
  res@cnLineLabelBackgroundColor = "transparent"
  res@cnMissingValFillColor        = "gray30"
  res@cnMissingValPerimOn          = True

  res@cnMissingValFillPattern      = -1;                set the missing value fill pattern to 5
  res@cnMissingValFillScaleF       = 0.9;              increase the density of the fill pattern (default   = 1.0)
;    resource@cnMissingValPerimOn  = True;             already turned on above
  res@cnMissingValPerimColor       = "black";          change the missing value perimeter to black
  res@cnMissingValPerimDashPattern = 1;           set the dash pattern of the missing value perimeter to 1
  res@cnMissingValPerimThicknessF  = 3.0;         increase the thickness of the missing value perimeter 3X

  res@gsnRightString	       = "[Sv]"
  res@gsnLeftString	       = "$basename"
  res@trYReverse = True	                 	; reverse the Y-axis
  res@tiXAxisString            = "latitude"
  res@tiYAxisString            = "depth in m"
  res@lbOrientation	       = "vertical" 
; res@lbLabelPosition          = "Left"  
; res@pmLabelBarOrthogonalPosF = 0.1           ; move label bar closer

; vertical levels with higher resolution near surface
  res@gsnYAxisIrregular2Linear = True          ; true: axis-scaling is linear in Y, false: following resolution

;*****************************************
;Experiment 014 (gi6)
;*****************************************


  delete(res@tmXBValues)
  delete(res@tmXBLabels)

  res@tmXBValues 	       = ispan(0,180,10)
  res@tmXBLabels 	       = ispan(-90,90,10)
; res@tiMainString	       = "MOC Global "  ; + infilename
; res@tiMainString	       = "MOC Global - $basename "  ; + infilename
  res@tiMainString	       = "MOC Global"

  plot = gsn_csm_contour(wks,moc_glo,res)

  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,120,10)
  res@tmXBLabels 	       = ispan(-30,90,10)
  res@tiMainString	       = "MOC Atlantic"
 plot = gsn_csm_contour(wks,moc_atl,res)

  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,110,10)
  res@tmXBLabels 	       = ispan(-50,70,10)
  res@tiMainString	       = "MOC Pac+Ind"
 plot = gsn_csm_contour(wks,moc_pac_ind,res)
end
EOF


ncl scr_plot_moc_my.ncl
rm scr_plot_moc_my.ncl
rm scr_moc_???.nc





