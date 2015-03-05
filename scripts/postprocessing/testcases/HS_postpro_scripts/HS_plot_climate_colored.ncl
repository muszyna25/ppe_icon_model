;
; Hui Wan ( MPI, Aug 2006 )
; Hui Wan ( MPI, Nov 2008 )
; Pilar Ripodas ( DWD, Oct 2010)
;---------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------
 begin

    PlotFormat = getenv("plot_file_format") 
    Model      = getenv("Model")
    DataPath   = getenv("DataPath")
    DataID     = getenv("DataID")
    DataGrid   = getenv("trunc")
    PlotPath   = getenv("PlotPath")
    ConfigStr  = getenv("ConfigStr")
    Resolution = getenv("Resolution")
    CellType   = getenv("cell_type")


    if (CellType .eq. 3) then
      cell_type=" triangles"
    end if

    if (CellType .eq. 6) then
      cell_type=" hexagons"
    end if
;------

 if (Model.eq."ICOHDC" .or. Model.eq."ICONAM" ) then

    suffix     = ".T"+DataGrid+".nc"

    varNameF   = (/"U.zmta","T.zmta","UpVp.zmta","TpVp.zmta","eddy-KE.zmta","TpTp.zmta"/) 
    varName    = (/"U",     "T",     "U",        "T",        "U",      "T"/)

    ;apzero     = 101325.  ; ECHAM5 default (ncar 2008: apzero= 1.e5)
    apzero     = 1.e5

;------
 else

 ;  suffix     = ".hyb.nc"

 ;  varNameF   = (/"ensavg.U","ensavg.T","ensavg.VpUp","ensavg.VpTp","ensavg.KE","ensavg.TT"/) 
 ;  varName    = (/"u",       "st",      "u",          "st",         "u",      "st"/)

 ;  apzero     = 101325.

 end if
;------

 LeftString = (/"(a) u-wind (ms~S~-1~N~)",      \
                "(b) temperature (K)", \
                "(c) eddy momentum flux (m~S~2~N~s~S~-2~N~)",    \
                "(d) eddy heat flux (Kms~S~-1~N~)",        \
                "(e) eddy kinetic energy (m~S~2~N~s~S~-2~N~)",   \
                "(f) temperature variance (K~S~2~N~)"/)

 nvar = dimsizes(LeftString)

 Min = (/ -12, 190., -90., -24,   0,  0/)
 Max = (/  32, 305.,  90.,  24, 480, 52/)
 Int = (/   4,   5,  15.,   3,  40,  4/)

 colorStart = (/  3,   3, 102, 102,   5,   5/)
 colorEnd   = (/101, 101, 201, 201, 101, 101/)

 nlev4plot = 101

;---------------------------------------------------------
; make plots
;---------------------------------------------------------

    wks = gsn_open_wks( PlotFormat,PlotPath+DataID+"_T"+DataGrid+"_climate")
    gsn_merge_colormaps( wks,"WhBlGrYeRe","testcmapshort")

    contour = new( nvar,graphic )

     ResC = True

     ResC@gsnDraw  = False
     ResC@gsnFrame = False

     if (Model .eq."ICOHDC") then
      ResC@tiYAxisString          = "Hybrid coordinate ~F8~h~F~" 
     end if
     if (Model .eq."ICONAM") then
      ResC@tiYAxisString          = "Height Km" 
     end if     

     FontHeight = 0.02

     ResC@tiXAxisFontHeightF        = FontHeight
     ResC@tiYAxisFontHeightF        = FontHeight
     ResC@tmXBLabelFontHeightF      = FontHeight
     ResC@tmYLLabelFontHeightF      = FontHeight
     ResC@gsnLeftStringFontHeightF  = FontHeight + 0.003
    ;ResC@gsnCenterStringFontHeightF= FontHeight + 0.004
     ResC@gsnRightStringFontHeightF = FontHeight + 0.002

     ResC@lbLabelBarOn             = True
     ResC@pmLabelBarHeightF        = 0.1
     ResC@pmLabelBarWidthF         = 0.58
     ResC@lbLabelAutoStride        = True
     ResC@pmLabelBarOrthogonalPosF = 0.01
     ResC@lbLabelFontHeightF       = FontHeight 

     ;ResC@trYLog        = False
     ResC@tmYRMode      = "Automatic"
     if (Model .eq."ICOHDC") then     
      ResC@trYReverse    = True
      ResC@trYMaxF       = 1.
      ResC@trYMinF       = 0.
     end if
     if (Model .eq."ICONAM") then 
      ;ResC@trGridType="LogLin"
      ;ResC@trYAxisType="LogAxis"
      ;ResC@trYLog        = True    
      ResC@trYReverse    = False
      ResC@trYMaxF       = 22.
      ResC@trYMinF       = 0.
      ;ResC@trYLog=True
      ;ResC@tmYLMinorPerMajor=8
     end if 

     ResC@cnInfoLabelOn                = False
   ; ResC@cnLineThicknessF             = 1.2
     ResC@cnLinesOn            = True
     ResC@gsnContourZeroLineThicknessF = 1.5
     ResC@gsnContourNegLineDashPattern = 2

     ResC@cnFillOn                     = True
     ResC@gsnSpreadColors              = True
    ;ResC@cnLineLabelsOn  = True
    ;ResC@cnLineLabelPlacementMode     = "Computed"
    ;ResC@cnLineLabelInterval          = 4
    ;ResC@cnLineLabelAngleF            = 0.0
    ;ResC@cnLineLabelFontHeightF       = 0.016
    ;ResC@cnLineLabelFontColor         = "black"
    ;ResC@cnLabelMasking               = True
    ;ResC@cnLineLabelBackgroundColor   = "transparent"
    ;ResC@cnLineLabelFormat          = "0@*+^.2g"



;---------------------------------------------------------
; data reading and plotting
;---------------------------------------------------------

 do ivar=0,nvar-1
    File = addfile( DataPath + DataID +"." + varNameF(ivar) +suffix, "r" )
    var = (/File->$varName(ivar)$(0,:,:,0)/)  ;(time,lev,cell,lon)

    if (ivar.eq.0) then
      if (Model .eq. "ICOHDC") then
       vertc = (File->hyam)/apzero + File->hybm ; eta as vertical coordinate
       vertc_t = 0.
       vertc_sfc = 1.
      end if
      if (Model .eq. "ICONAM") then
       height = (File->hyam)      ; orography ist 0. in this test case
       vertc =  height/1000.      ; height in km as vertical coordinate
       vertc_t = 22.
       vertc_sfc = 0.
      end if
    end if

    var!0   = "vertc"
    var&vertc = vertc
    var!1   = "lat"
    var&lat = File->lat

    vertc4plot     = fspan(vertc_t,vertc_sfc,nlev4plot)
    vertc4plot!0   = "vertc"
    vertc4plot&vertc = vertc4plot

    tmp      = linint1_Wrap( vertc, var(lat|:,vertc|:), False, vertc4plot, 0 )
    var4plot = tmp(vertc|:,lat|:)

     ResC@cnLevelSelectionMode  = "ManualLevels"
     ResC@cnMinLevelValF        =  Min(ivar)
     ResC@cnMaxLevelValF        =  Max(ivar)
     ResC@cnLevelSpacingF       =  Int(ivar)

     ResC@gsnLeftString  = LeftString(ivar)
     ResC@gsnRightString = Resolution ;ConfigStr

     ResC@gsnSpreadColorStart = colorStart(ivar)
     ResC@gsnSpreadColorEnd   = colorEnd(ivar)

     contour(ivar) = gsn_csm_contour( wks,var4plot,ResC )
  end do

     ResP                            = True
     ResP@gsnMaximize                = True
   ; ResP@gsnPaperMargin             = 0.2
     ResP@gsnPanelXWhiteSpacePercent = 5
     ResP@gsnPanelYWhiteSpacePercent = 5
     ResP@txString = DataID+cell_type
     ResP@txFontHeightF  = 0.01

     gsn_panel (wks, contour, (/3,2/),ResP)

 end
