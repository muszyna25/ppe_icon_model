#!/bin/bash
#

# plot T,S and potential density variation to initial values from a list of ICON input files

#  updates Ralf
#   - ncl (icon_plot) with dozens of errors on blizzard, but okay on workstation
#   - cdo -r cat necessary - done

#  updates Stephan:
#   - use all timesteps available -> removes bias at decades/eof-file
#   - use new variable names t_acc, s_acc, rhopot_acc of 'nml'-output
#   - an init file is generated using 10 days output (2001-01-11), since first timestep contains zero

#  TBD:
#

# ==============================================================================
set -ex
#  Tools and their paths
     CDO="cdo"                                   #  climate data operators
ICONPLOT=/pool/data/ICON/tools/icon_plot.ncl     #  NCL-script
 ICONLIB=/pool/data/ICON/tools

# ==============================================================================
# ==============================================================================
# we supppose, that all files belong to the same experiment
         expIdent='xmpiom.r12088.norunoff'                                           # experiment identifier
         expIdent='xmpiom.r12690.upw.d400'                                           # experiment identifier
         expIdent='xmpiom.r13694'                                                    # experiment identifier
          expPath='/scratch/mpi/mh0287/users/m211032/Icon/icon-dev.new/experiments'  # experiment path
     fileListPath="$expPath/$expIdent/Output"                                        # output data path
     fileListPath="$expPath/$expIdent"                                               # output data path
  fileListPattern="${expIdent}_iconR2B04-ocean_etopo40_planet_000[1-5].nc"           # 'nml' naming convention
  fileListPattern="${expIdent}_R2B04_oce_DOM01_ML_000[1-5].nc"                       # 'vlist' naming convention
      outputIdent='r12088.relice'                                                    # output file name appendix
      outputIdent='r13694'                                                           # output file name appendix
          Tempvar='T'                                                                # temperature variable name
           Salvar='S'                                                                # salinity variable name
           Rhovar='rhopot'                                                           # density variable name
          Tempvar='t_acc'                                                            # temperature variable name
           Salvar='s_acc'                                                            # salinity variable name
           Rhovar='rhopot_acc'                                                       # density variable name
   outputDataFile="ano.TSrho.$outputIdent.nc"                                        # output data file name
       PlotScript='ncl'               #  ncl-script, see below - not yet
       PlotScript='icon'              #  script icon_plot_ncl using ncl, see above
       PlotScript='none'              #  no plots
# ==============================================================================
# these variable names must not necessarily be changed
      MaskVarName='wet_c'
          ScrDatN='tsrho'
          ScrDatN='scr_tsrho'
# ==============================================================================
#declare a fileListArray
fileList=$(ls $fileListPath/$fileListPattern)
i=0
for file in $fileList; do
  fileListArray[$i]=$file
  i=$((i+1))
done
numberOfFiles=$i
echo "$numberOfFiles number of files"
echo $fileList
# ==============================================================================
# handing of the initial values
# get the initial values before any averaging is done
#initFile='init-phc.r2b4.20L.ts_acc.nc'     #  orig PHC T/S data
# include rhopot in PHC data:
# cdo chname,rhopoto,rhopot_acc -rhopot -adisit -chname,t_acc,tho,s_acc,sao init-phc.r2b4.20L.ts_acc.nc xx
#  - now use first data set of output file for subtracting init!
#    otherwise sequence, names, dimensions etc. are different when using PHC data
#  - caution: first set of accumulated values is zero
#initFile='init-phc.rho.r2b4.20L.ts_acc.nc'  #  contains names as in ocean output: t_acc, s_acc, rhopot_acc
initFile=fldmean_init.$outputIdent.nc

# ==============================================================================
# using a mask file - contains no missing values
maskFile='wetc.r2b4.r11xxx.nc'    # easy to generate - cdo selvar,wet_c
maskFile="wetc.$expIdent.nc"      # current mask created automatically
# ==============================================================================
# Loop over all files in serial
for file in $fileList; do
  fnum=${file##*_}  #  '_' must be last part before "0001.nc" of ocean output filename
  pattern_num=${outputIdent}_$fnum
  if [[ -f fldmean_${pattern_num} ]]; then
    echo  "fldmean_${pattern_num} already exists!"
  else
    if [[ -f ${ScrDatN}_${pattern_num} ]]; then
      echo  "${ScrDatN}_${pattern_num} already exists!"
    else
      if [[ ! -f "$maskFile" ]]; then 
        # create maskFile
        $CDO selvar,$MaskVarName $file $maskFile
      fi
      # select variables, mask out land points - takes some time
      $CDO -div -selname,$Tempvar,$Salvar,$Rhovar $file $maskFile ${ScrDatN}_${pattern_num}
    fi
    # compute the fldmean using masked values only
    $CDO fldmean ${ScrDatN}_${pattern_num} fldmean_${pattern_num}

    # special treatment of first file "*0001.nc": contains zero!
    #  - for _acc accumulated variables only
    if [[ $Tempvar == "t_acc" ]]; then 
      if [[ $fnum == "0001.nc" ]]; then 
        # no deldate available?
        # seldate,2001-01-11,2010-12-30 for 10-days output only!
        $CDO seldate,2001-01-11,2010-12-30 fldmean_${pattern_num} fldmean.noinit.nc
        mv fldmean.noinit.nc fldmean_${pattern_num}
      fi  # file no=1
    fi  # tempvar=t_acc
  fi  # create fldmean_${pattern_num}
  if [[ $fnum == "0001.nc" ]]; then 
    # now use first timestep for hovmoeller anomaly plot:
    $CDO seltimestep,1 fldmean_${pattern_num} $initFile
  fi
done

# ==============================================================================
# Cat the files together - results would be appended to last file, delete it
[[ -f fldmean_$outputIdent.nc ]] && rm fldmean_$outputIdent.nc
$CDO -r cat fldmean_${outputIdent}_* fldmean_$outputIdent.nc

# ==============================================================================
# Subtract initial values
[[ -f $outputDataFile ]] || $CDO yearmean -sub fldmean_$outputIdent.nc $initFile $outputDataFile

# ==============================================================================
# Plot a hovmoeller type graph using icon_plot.ncl
#  reports lots of errors, difficult to put on one page (psnup -3 if.ps of.ps)
if [[ $PlotScript == "icon" ]]; then 
  for varname in $Tempvar $Salvar $Rhovar; do 
    nclsh $ICONPLOT  -altLibDir=$ICONLIB -varName=$varname -iFile=$outputDataFile \
    -oFile=hov.ano.${outputIdent}_$varname -oType=ps -isIcon -hov=true
  done
fi

# ==============================================================================
# Plot a hovmoeller type graph using ncl-script
if [[ $PlotScript == "ncl" ]]; then 

cat >scr_plot_tsrho.ncl << EOF
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************

begin

;================================================
; get data
;================================================
  f     = addfile ("$outputDataFile", "r")
  temp  = f->$Tempvar(:,:,:,:)

; t     = f->$Tempvar(depth|time|lat|lon)   ; error ??
; s     = f->$Salvar(:,:,:,:)
; rho   = f->$Rhovar(:,:,:,:)

; t     = temp(depth|:,time|0|0)   ; wrong sequence
;moc_atl     = f->var778(0,:,60:179)  ; 120 values - 30S to 90N
; t_    = t(depth_2|:,time|:,lat|:,lon|:)
; t     = f->$Tempvar(:,:,0,0)

  t     = temp(time|:,depth|:,lat|:,lon|:)


;************************************************************************
;  more TBD:
;   - ???
;************************************************************************

  wks = gsn_open_wks("ps","hov.ano.${outputIdent}")
; wks = gsn_open_wks("pdf","hov.ano.${outputIdent}")

  aplot = new(3,graphic);  resources for 3 plots on one page

; gsn_define_colormap(wks,"BlWhRe")            ; colormap
; gsn_define_colormap(wks,"BlGrYeOrReVi200")   ; colormap
; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
; gsn_define_colormap(wks,"temp_diff_18lev")
; gsn_define_colormap(wks,"temp_19lev")
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
; gsn_define_colormap(wks,"BlueWhiteOrangeRed")
; gsn_define_colormap(wks,"testcmap")
  gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

  res                      = True
  res@gsnDraw              = False             ; False: don't draw on own page
  res@gsnFrame             = False             ; False: don't advance frame - otherwise empty plots?

; res@cnFillMode           = "RasterFill"
  res@gsnPaperOrientation  = "Portrait"
; res@gsnPaperOrientation  = "Landscape"
; res@gsnMaximize          = True              ; no effect

; res@gsnYAxisIrregular2Linear = True          ; True: axis-scaling is linear in Y, False: following resolution

  res@vpWidthF             = 0.85              ; change aspect ratio of plot
  res@vpHeightF            = 0.32              ; change aspect ratio of plot

  res@cnFillOn             = True              ; color for contours
; res@cnLinesOn            = False             ; contour lines
; res@lbLabelBarOn         = False             ; color label bar
  res@cnInfoLabelOn        = False             ; contour info label
  res@gsnSpreadColors      = True              ; subset of the whole colormap
; res@gsnSpreadColorStart  =  25               ; start color
; res@gsnSpreadColorEnd    = -25               ; end color
; res@cnFillMode           = "CellFill"        ; filled cells


; choice: explicit levels
; res@cnLevelSelectionMode  = "ExplicitLevels"  ; set explicit contour levels
; res@cnLevels              = (/ -100, -20, -10, -5,0,5,10,20,100 /)

  res@gsnContourZeroLineThicknessF   = 2.2      ; factor for thickness of zero contour
  res@gsnContourLineThicknessesScale = 1.5      ; factor for other lines

; res@cnLineLabelsOn        = False             ; contour line labels
; res@cnLineLabelFontHeightF= 0.012             ; contour line label height
; res@cnLineLabelInterval   = 2                 ; line label interval (default=2)
; res@cnLabelMasking        = True		; masking contour lines below line labels
; res@cnLineLabelBackgroundColor = "white"

  res@lbLabelFontHeightF    = 0.014             ; label bar font
; res@lbLabelStride         = 2                 ; labelling interval at color bar
; res@lbLabelPosition       = "Left"  
; res@pmLabelBarOrthogonalPosF = 0.1            ; move label bar closer

; res@cnMissingValFillColor        = "gray30"    ; not defined?
; res@cnMissingValPerimOn          = True
; res@cnMissingValFillPattern      = -1          ; set the missing value fill pattern
; res@cnMissingValFillScaleF       = 0.9         ; increase the density of the fill pattern (default   = 1.0)
; res@cnMissingValPerimColor       = "black"     ; change the missing value perimeter to black
; res@cnMissingValPerimDashPattern = 1           ; set the dash pattern of the missing value perimeter to 1
; res@cnMissingValPerimThicknessF  = 2.0         ; factor for thickness of the missing value perimeter 3X

  res@gsnRightString	        = "hov.ano.${outputIdent}"  ;  output filename
  res@trYReverse                = True            ; reverse the Y-axis
                            
  res@tiXAxisOn                 =  True        ;  false: x-axis title removed - cannot be redrawn?
; res@tiXAxisString             = "Time [years]"
  res@tiXAxisString             = " "          ;  set to empty, used only at bottom
  res@tiXAxisFontHeightF        = 0.025

  res@tiYAxisString             = "depth [m]"
  res@tiYAxisFontHeightF        = 0.025

; resources used by KMF - to be checked
; res@tmXTOn                    =  True        ;  false: turns off top ticks
; res@tmXTLabelsOn              =  True        ;  false: turns off top labels
; res@tmXBOn                    =  False       ;  false: x-axis numbers and ticks removed
  res@tmXBMode 		        = "Explicit"   ;  switch from dimension to latitude
; res@tmXBLabelsOn              =  False       ;  false: x-axis labels removed
  res@tmXBLabelFontHeightF      = 0.020        ;  x-axis labels font height
  res@tmYLLabelFontHeightF      = 0.018        ;  y-axis labels font height

; *****************************************
; Plotting the 3 panels
; *****************************************

; choice: range of levels
  res@gsnLeftString	    = "$Tempvar [K]"
  res@cnLevelSelectionMode  = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF        = -3                ; set min contour level
  res@cnMaxLevelValF        = 3                 ; set max contour level
  res@cnLevelSpacingF       = 0.3               ; set contour spacing

; drawNDCGrid(wks)  ;  - verkleinert
  plot     = gsn_csm_contour(wks,t(:,:,0,0),res)
  aplot(0) = plot   ;  creates the uppermost plot

; ==========  ========== ==========

  res@gsnLeftString	    = "$Salvar [psu]"
  res@cnMinLevelValF        = -0.3              ; set min contour level
  res@cnMaxLevelValF        = 0.3               ; set max contour level
  res@cnLevelSpacingF       = 0.03              ; set contour spacing

; plot     = gsn_csm_contour(wks,s,res)
; aplot(1) = plot

; ==========  ========== ==========

  res@gsnLeftString	    = "$Rhovar [kg/m3]"
  res@cnMinLevelValF        = -0.5              ; set min contour level
  res@cnMaxLevelValF        = 0.5               ; set max contour level
  res@cnLevelSpacingF       = 0.05              ; set contour spacing

; res@tiXAxisOn                = True         ;  false: x-axis title removed
  res@tiXAxisOffsetYF          =  0.010 
  res@tiXAxisString            = "time [years]"
; plot     = gsn_csm_contour(wks,rho,res)
; aplot(2) = plot

; *****************************************
;  Resources for 3 panels
; *****************************************

; Attention
;  - the maximize function (and other carackteristics of the whole page using gsn_panel) need own resources (resP)
  resP                           = True
  resP@gsnPaperOrientation       = "Portrait"
; resP@gsnPaperOrientation       = "Landscape"
  resP@gsnMaximize               =  True
  resP@txString                  = "ICON: Global Mean Evolution"
  resP@txFontHeightF             = 0.020
; resP@tiMainString	         = "ICON: Meridional Overturning Stream Function"
; resP@tiMainOffsetYF            =  -0.3       ;  wo ist der title
; resP@tiMainOn                  = True        ;  False: Main title removed

; resP@gsnPanelYWhiteSpacePercent =  0.1       ;  <-- KMF
; resP@gsnFrame                  =  False
; resP@gsnPaperHeight            =  11.69      ;  <-- KMF
; resP@gsnPaperWidth             =   8.27      ;  <-- KMF

  resP@gsnPanelBottom            =  0.02       ;  <-- KMF
  resP@gsnPanelTop               =  0.95       ;  <-- KMF
; resP@tiXAxisOn                 =  False      ;  <-- KMF
; resP@vpXF                      =  0.05
; resP@vpYF                      =  0.95

; drawNDCGrid(wks)                       ;  - verkleinert Plot?
  gsn_panel(wks,(/aplot/),(/3,1/),resP)  ; now draw as one plot


end
EOF


  ncl scr_plot_tsrho.ncl
# rm scr_plot_tsrho.ncl


fi


