;------------------------------------
; Hui Wan (MPI-Met, 2010-07-30)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes contour plots of surface pressure and 
; some precipitation related quantities for the 
; Jablonowski-Williamson wave state test with ECHAM physics.
;---------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;---------------------------------------------------------------
 begin

  PlotFormat = getenv("plot_file_format") 
  Model      = getenv("Model")
  DataPath   = getenv("DataPath")
  PlotPath   = getenv("PlotPath")
  Resolution = getenv("Resolution") 
  ConfigStr  = getenv("ConfigStr")
  ExpName    = getenv("ExpName")
  plotBaseName = getenv("plotBaseName")
  CellType   = getenv("cell_type")
  OutFrq     = stringtoint(getenv("output_frequency"))
  dayStart   = stringtoint(getenv("post_start"))
  dayEnd     = stringtoint(getenv("post_end"))
  multiPanel = stringtoint(getenv("multi_panel"))

  rad2deg = 45./atan(1.)   ; radians to degrees

  colormap   = "WhViBlGrYeOrReWh"
  colorStart = 2
  colorEnd   = 100

;---------------------------------------------------------------
; time steps to be plotted
;---------------------------------------------------------------

  iday  = ispan(dayStart,dayEnd,1)  ; days to plot
  ndays = dimsizes(iday)

  tid = iday
  if (Model.eq."ECHAM") then
     tid = 0
  else
     tid = OutFrq - 1
     if (dayStart.eq.1) then
        tid(0) = OutFrq   ; 1st file contains also the initial conditions
     end if
  end if

;---------------------------------------------------------------
; data and plot files; variables to plot
;---------------------------------------------------------------

  if (Model.eq."ECHAM") then
     DataFileName = DataPath + ExpName +"_"+sprinti("%4.4i",iday+1)+".nc"
  else
     DataFileName = DataPath + ExpName +"_"+sprinti("%4.4i",iday)+".nc"
  end if
  PlotFileName = PlotPath + plotBaseName +"_"+sprinti("%4.4i",iday)+"_precip"

  if (Model.eq."ICOHAM") then
     prName = (/"APRL","RSFL","SSFL","APRC","RSFC","SSFC"/)
     psName = "PS"
  else
     prName = (/"aprl","jrsfl","jssfl","aprc","jrsfc","jssfc"/)
     psName = "aps"
  end if

  prMin = (/ 1e-5,  2e-5,  2e-5,  1e-5,  2e-5,  2e-5/)
  prMax = (/15e-5, 30e-5, 30e-5, 15e-5, 30e-5, 30e-5/)
  prInt = (/ 1e-5,  2e-5,  2e-5,  1e-5,  2e-5,  2e-5/)

  VarLongName = prName
  unit        = prName
  unit        = "Kg m~S~-2~N~ s~S~-1~N~" 

  nvar = dimsizes(prName)
  plot = new(nvar,graphic)
  dum  = plot

;---------------------------------------------------------------
; lat./lon. of cell center (ICON only)
;---------------------------------------------------------------

  if (Model.eq."ICOHAM") then

     if (CellType.eq."3") then
        latName    = "clat"
        lonName    = "clon"
     else
     if (CellType.eq."6") then
        latName    = "vlat"
        lonName    = "vlon"
     else
        print("Wrong cell type. Should be 3 or 6")
     end if
     end if

     File = addfile( DataFileName(0), "r" )
     x = File->$lonName$ *rad2deg
     y = File->$latName$ *rad2deg
  end if

;---------------------------------------------------------------
; common plot settings
;---------------------------------------------------------------

  ResC                = True
  ResC@gsnDraw        = False
  ResC@gsnFrame       = False
  ResC@cnLineLabelsOn = False                                   

  FontHeight = 0.011                                            
                                                                
  ResC@tiXAxisFontHeightF    = FontHeight                       
  ResC@tiYAxisFontHeightF    = FontHeight                       
  ResC@tmXBLabelFontHeightF  = FontHeight                       
  ResC@tmYLLabelFontHeightF  = FontHeight                       
  ResC@gsnStringFontHeightF  = FontHeight +0.002                

  if (Model.eq."ICOHAM") then
     ResC@sfXArray = x
     ResC@sfYArray = y
  end if

  ResL = ResC

; shadings vs contours                                          

  ResC@cnFillOn       = True                                    
  ResC@cnLinesOn      = False                                   
  ResC@cnInfoLabelOn  = False                                   
                                                                
  ResL@cnFillOn         = False                                  
  ResL@cnLinesOn        = True
  ResL@cnLineThicknessF = 0.6     
  ResL@cnInfoLabelOn             = True                                   
  ResL@cnInfoLabelFontHeightF    = FontHeight*0.7
  ResL@cnInfoLabelOrthogonalPosF = -0.55
  ResL@cnInfoLabelParallelPosF   = 0.99
  ResL@cnInfoLabelPerimOn        = False

; label bar (only) for the color plot                             
                                                                
  ResC@lbLabelBarOn             = True              
  ResC@pmLabelBarHeightF        = 0.06                          
  ResC@pmLabelBarWidthF         = 0.75                           
  ResC@lbLabelAutoStride        = True                          
  ResC@pmLabelBarOrthogonalPosF = 0.15                      
  ResC@lbLabelFontHeightF       = FontHeight                    
                                                                
  ResL@lbLabelBarOn             = False              

; contour levels                                                
                                                                
  ResC@cnLevelSelectionMode = "ManualLevels"                    
  ResC@gsnSpreadColors      = True                              
  ResC@gsnSpreadColorStart  = colorStart                        
  ResC@gsnSpreadColorEnd    = colorEnd                          
                                                                
  ResL@cnLevelSelectionMode = "ManualLevels"                   

; the color shading plot is a "map plot"                        
                                                                
  ResC@mpProjection           = "CylindricalEquidistant"        
  ResC@mpGeophysicalLineColor = "transparent"                   
  ResC@mpFillOn               = False                           
  ResC@mpLimitMode            = "LatLon"                        
  ResC@mpMinLatF              = 0.                              
  ResC@mpMaxLatF              = 90.                             
  ResC@mpMinLonF              = 0.                              
  ResC@mpMaxLonF              = 360.                            
  ResC@mpCenterLonF           = 180.                            

  ResC@gsnMajorLatSpacing = 30.
  ResC@gsnMajorLonSpacing = 60.
  ResC@gsnMinorLatSpacing = 10.
  ResC@gsnMinorLonSpacing = 10.

  ResL@gsnLeftString   = "" 
  ResL@gsnCenterString = "" 
  ResL@gsnRightString  = "" 

;---------------------------------------------------------------
; read data and make plots
;---------------------------------------------------------------
  do it = 0,ndays-1

     File = addfile( DataFileName(it), "r" )

     if (Model.eq."ICOHAM") then
        ps = File->$psName$(tid(it),:)
     else ;ECHAM
        ps = File->$psName$(tid(it),:,:)
     end if
     ps = ps*0.01

     if (multiPanel.eq.1) then
        wks = gsn_open_wks(PlotFormat,PlotFileName(it))
        gsn_define_colormap(wks,colormap)
     end if

     ResC@gsnCenterString = "day "+iday(it) 
     ResC@gsnRightString  = ConfigStr

     if (iday(it).lt.7) then
       ;ResL@cnMinLevelValF       =  993.
       ;ResL@cnMaxLevelValF       = 1005.
        ResL@cnLevelSpacingF      =    1. 
     else
        ResL@cnMinLevelValF       = floattointeger(min(ps)/10)*10
        ResL@cnMaxLevelValF       = floattointeger(max(ps)/10)*10
        ResL@cnLevelSpacingF      =   10.
     end if

     do iv = 0,nvar-1

        if (multiPanel.ne.1) then  ; one panel per file
           wks = gsn_open_wks(PlotFormat,PlotFileName(it)+"_"+prName(iv))
           gsn_define_colormap(wks,colormap)
        end if

        if (Model.eq."ICOHAM") then
           pr = File->$prName(iv)$(tid(it),:)
           if (iv.eq.0.or.iv.eq.3) then
             pr = pr*OutFrq/86400.
           end if
        else ;ECHAM
           pr = File->$prName(iv)$(tid(it),:,:)
        end if

        ResC@gsnLeftString   = "ps (hPa, lines) and +"+ VarLongName(iv)+" ("+unit(iv)+", colors)"
        ResC@cnMinLevelValF  = prMin(iv)
        ResC@cnMaxLevelValF  = prMax(iv)
        ResC@cnLevelSpacingF = prInt(iv)

        plot(iv) = gsn_csm_contour_map(wks,pr,ResC)
        dum (iv) = gsn_csm_contour    (wks,ps,ResL)
        overlay( plot(iv),dum(iv) )

        if (multiPanel.ne.1) then  ; one panel per file
           draw(plot(iv))
           frame(wks)
           delete(wks)
        end if
     end do

     if (multiPanel.eq.1) then  ; one panel per file
        ResP = True
        ResP@gsnMaximize = True
        ResP@gsnPanelYWhiteSpacePercent = 5
        ResP@txString = "Precipitation (color) and surface pressure (contour)"
        ResP@txFontHeightF = FontHeight
        gsn_panel(wks,plot,(/nvar,1/),ResP)
        delete(wks)
     end if

     print(" day "+iday(it)+" done")
  end do
end

