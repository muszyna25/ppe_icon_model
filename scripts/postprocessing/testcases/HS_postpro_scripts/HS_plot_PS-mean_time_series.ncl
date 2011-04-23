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
  Exp        = getenv("EXP")

  title = "Held-Suarez Test"
  DataFile = DataPath+DataID+".nc" 
  PlotFile = PlotPath+DataID

  if (CellType .eq. 3) then
    cell_type="  triangles"
  end if

  if (CellType .eq. 6) then
    cell_type="  hexagons"
  end if
;----------------------------------------------


  File = addfile( DataFile,"r" )

;--------------------------------------------------------------- 
; plot settings

     ResXY = True

     ResXY@gsnMaximize    = True
     ResXY@gsnPaperOrientation="portrait"

     ResXY@tiYAxisString          = "Surface Pressure Globa Mean (hPa)"
     ResXY@tiYAxisFontHeightF     = 0.016
     ;ResXY@trYMinF                = 950.
     ;ResXY@trYMaxF                = 1005.

     ResXY@tiXAxisString          = "Time ( days )"
     ResXY@tiXAxisFontHeightF     = 0.016



     wks = gsn_open_wks( PlotFormat,PlotFile )

     var = File->$VarName$(:,0,0)  ;(time,lat,lon)
     var = var/100.

     ntime = dimsizes(var&time)
     
     time = fspan(0,ntime-1,ntime)

     ResXY@gsnLeftString   = "" ;ConfigStr
     ResXY@gsnCenterString = ""
     ResXY@gsnRightString  = "" ;"Evolution of Surface Pressure global mean" 
     ResXY@tiMainString = Exp+"  "+Resolution+cell_type
     ResXY@tiMainFontHeightF  = 0.015


     contour = gsn_csm_xy( wks,time,var,ResXY )

end
