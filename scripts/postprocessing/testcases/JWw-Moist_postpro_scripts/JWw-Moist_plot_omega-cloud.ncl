;------------------------------------
; Hui Wan (MPI-Met, 2008-07-13)
; Hui Wan (MPI-Met, 2009-04-28: Modification for the hex. model)
; Hui Wan (MPI-Met, 2010-07-30)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes contour plots of 850 hPa temperature 
; for the Jablonowski-Williamson steady state test.
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
     DataFileName = DataPath + ExpName +"_"+sprinti("%4.4i",iday+1)+"_cloud_p.nc"
  else
     DataFileName = DataPath + ExpName +"_"+sprinti("%4.4i",iday)+"_cloud_p.nc"
  end if
  PlotFileName = PlotPath + ExpName +"_"+sprinti("%4.4i",iday)+"_cloud"

  if (Model.eq."ICOHAM") then
     cldName = (/"ACLC","Qw","Qi"/)
     omgName = "OMEGA_PHY"
  else 
  if (Model.eq."ECHAM") then
     cldName = (/"aclc","xl","xi"/)
     omgName = "var135"
  else
     print(" Wrong model name! Abort.")
     exit
  end if
  end if

  VarLongName = cldName
  unit        = cldName
  unit        = "" 

  nvar = dimsizes(cldName)
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
  else

     File = addfile( DataFileName(0), "r" )
  end if

;---------------------------------------------------------------
; pressure levels
;---------------------------------------------------------------

  lev = File->lev
  if (lev@units.eq."Pa") then
     lev = lev*0.01
  end if
  levunit = "hPa"
  nlev = dimsizes(lev)

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
  ResL@gsnContourZeroLineThicknessF = 1.
  ResL@gsnContourNegLineDashPattern = 2
  ResL@cnLineDashSegLenF            = 0.05


; label bar (only) for the color plot                             
                                                                
  ResC@lbLabelBarOn             = True              
  ResC@pmLabelBarHeightF        = 0.06                          
  ResC@pmLabelBarWidthF         = 0.75                           
  ResC@lbLabelAutoStride        = True                          
  ResC@pmLabelBarOrthogonalPosF = 0.15                      
  ResC@lbLabelFontHeightF       = FontHeight                    
                                                                
  ResL@lbLabelBarOn             = False              

; contour levels                                                
                                                                
 ;ResC@cnLevelSelectionMode = "ManualLevels"                    
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
  do jk = 0,nlev-1

     File = addfile( DataFileName(it), "r" )

     if (Model.eq."ICOHAM") then
        omg = File->$omgName$(tid(it),jk,:)
     else ;ECHAM
        omg = File->$omgName$(tid(it),jk,:,:)
     end if

     if (multiPanel.eq.1) then
        FNameTmp = PlotFileName(it) + lev(jk) + levunit
        wks = gsn_open_wks(PlotFormat,FNameTmp)
        gsn_define_colormap(wks,colormap)
     end if

     ResC@gsnCenterString = "day "+iday(it) 
     ResC@gsnRightString  = ConfigStr

     if (iday(it).lt.7) then
        ResL@cnMinLevelValF  = floattointeger(min(omg)*100)/100.
        ResL@cnMaxLevelValF  = floattointeger(max(omg)*100)/100.
        ResL@cnLevelSpacingF = 0.1
     else
        ResL@cnMinLevelValF  = floattointeger(min(omg)*10)/10.
        ResL@cnMaxLevelValF  = floattointeger(max(omg)*10)/10.
        ResL@cnLevelSpacingF = 0.2
     end if

     do iv = 0,nvar-1

        if (multiPanel.ne.1) then  ; one panel per file
          FNameTmp = PlotFileName(it) + lev(jk) + levunit +"_"+cldName(iv)
          wks = gsn_open_wks(PlotFormat,FNameTmp)
          gsn_define_colormap(wks,colormap)
        end if

        if (Model.eq."ICOHAM") then
           cld = File->$cldName(iv)$(tid(it),jk,:)
        else ;ECHAM
           cld = File->$cldName(iv)$(tid(it),jk,:,:)
        end if

        ResC@gsnLeftString   = lev(jk)+" "+levunit+" "+VarLongName(iv)+" and omega (Pa/s)"

        plot(iv) = gsn_csm_contour_map(wks,cld,ResC)
        dum (iv) = gsn_csm_contour    (wks,omg,ResL)
        overlay( plot(iv),dum(iv) )
        if (multiPanel.ne.1) then
           draw(plot(iv))
           frame(wks)
           delete(wks)
        end if 
     end do

     if (multiPanel.eq.1) then
        ResP = True
        ResP@gsnMaximize = True
        ResP@gsnPanelYWhiteSpacePercent = 5
        ResP@txString = "Cloud (color) and vertical velocity (contour)"
        ResP@txFontHeightF = FontHeight
        gsn_panel(wks,plot,(/nvar,1/),ResP)
        delete(wks)
     end if

     print(" day "+iday(it)+" level "+lev(jk)+" done")
  end do
  end do
end

