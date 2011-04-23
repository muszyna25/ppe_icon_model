;
; Hui Wan (MPI, 2005.07.27)
;-----------------------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------------------------------
 begin
  PlotFormat = getenv("plot_file_format")
  Model      = getenv("Model")
  DataPath   = getenv("DataPath")
  DataID     = getenv("DataID")
  PlotPath   = getenv("PlotPath")
  ConfigStr  = getenv("ConfigStr")
  VarName    = getenv("VarName")
  CellType   = getenv("cell_type")
  Resolution = getenv("Resolution")

  title = "Held-Suarez Test"
  DataFile = DataPath+DataID+".nc" 
  PlotFile = PlotPath+DataID

  if (CellType .eq. 3) then
    cell_type="tri"
  end if

  if (CellType .eq. 6) then
    cell_type="hex"
  end if
;----------------------------------------------


  File = addfile( DataFile,"r" )
  if ( Model .eq. "ICOHDC") then
     lev  = doubletointeger(File->lev)
     lev  = lev/100
  end if
  if ( Model .eq. "ICONAM") then
     lev  = File->lev
  end if
;--------------------------------------------------------------- 
; plot settings

     ResC = True

     ResC@gsnMaximize    = True
     ResC@vpHeightF      = 0.2
     ResC@vpWidthF       = 0.6
     ResC@vpXF           = 0.03
     ResC@vpYF           = 0.97

  ;  ResC@tiYAxisString          = "Latitude"
  ;  ResC@tiYAxisFontHeightF     = 0.016
     ResC@tmYLLabelFontHeightF   = 0.014
     ResC@tmYLMajorLengthF       = 0.01
     ResC@tmYLMinorLengthF       = 0.005
     ResC@tmYRMajorLengthF       = 0.01
     ResC@tmYRMinorLengthF       = 0.005
     ResC@tmYLMajorThicknessF    = 3. 
     ResC@tmYRMajorThicknessF    = 3.

     ResC@tiXAxisString          = "Time ( days )"
     ResC@tiXAxisFontHeightF     = 0.016

     ResC@tmXBLabelFontHeightF   = 0.014
     ResC@tmXBMajorLengthF       = 0.01
     ResC@tmXBMinorLengthF       = 0.005
     ResC@tmXBMajorThicknessF    = 3. 
     ResC@tmXTMajorThicknessF    = 3. 

     ResC@cnFillOn       = True
     ResC@cnLinesOn      = False
     ResC@cnInfoLabelOn  = False

     ResC@lbLabelBarOn             = True
     ResC@lbOrientation            = "vertical"
  ;  ResC@pmLabelBarOrthogonalPosF = -0.001
     ResC@pmLabelBarWidthF         = 0.07
     ResC@lbLabelFontHeightF       = .012
     ResC@lbLabelAutoStride        = True
     ResC@gsnSpreadColors          = True            ; use full range of colors
     ResC@gsnSpreadColorStart      = 2
     ResC@gsnSpreadColorEnd        = 176

     wks = gsn_open_wks( PlotFormat,PlotFile )
     gsn_define_colormap( wks,"BkBlAqGrYeOrReViWh200" )

     var = File->$VarName$(:,0,:,0)  ;(time,lev,lat,lon)

     ntime = dimsizes(var&time)
     delete(var&time)
     var&time = fspan(0,ntime-1,ntime)

     ;ResC@gsnLeftString   = "" ;ConfigStr
     ResC@gsnLeftString   = Resolution+cell_type
     ResC@gsnCenterString = ""
     if ( Model .eq. "ICOHDC") then
      ResC@gsnRightString  = "Evolution of zonal temperature variance (K~S~2~N~) at "+lev+" hPa"
     end if
     if ( Model .eq. "ICONAM") then
      ResC@gsnRightString  = "Evolution of zonal temperature variance (K~S~2~N~) at model level "+lev
     end if
     
     ResC@cnLevelSelectionMode  = "ManualLevels"
     ResC@cnMinLevelValF        =  0             
     ResC@cnMaxLevelValF        =  80 
     ResC@cnLevelSpacingF       =  5 

     contour = gsn_csm_lat_time( wks,var(lat|:,time|:),ResC )

end
