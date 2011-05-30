;---------------------------------------------------------------
; This script makes contour plots of ocean state
;---------------------------------------------------------------
;  - last changes by Stephan - 05/2011
;
;---------------------------------------------------------------
; Usage:
;   ncl icon_ocean.ncl 'iFile="path/file.nc"' 'varName="W"'    'mapLLC=(/35.0, -8/)' 'mapURC=(/55, 8/)'
;   ncl icon_ocean.ncl 'iFile="path/file.nc"' 'varName="ELEV"' 'oFile="test"' timeStep=1 levIndex=1

; Masking:
;   ncl icon_ocean.ncl 'iFile="oce_aqua_planet_O.R2B04L4_0001.nc"' 'varName="ELEV"' 'maskName="wet_c"' 
;   ncl icon_ocean.ncl 'iFile="oce.nc"' 'varName="W"' 'maskName="topo"' 'maskFile="icon_mask.nc"' 
;
; Vectorplot:
;   ncl icon_ocean.ncl 'iFile="oce.nc"' 'vecVars="u-veloc v-veloc"' 'oFile="test"'
;---------------------------------------------------------------
;
; Required Parameter:
; iFile         : input file
; oFile         : plot file
;
; Optional Parameter:
; oType         : output graphics format (ps, eps, png, ...)
; varName       : name of the variable to (contour) plot
; vecVars       : space separated string of the 2 vector components u and v to draw a vector plot
; resolution    : resolution string whish is used for remapping the icon data to regular grid
; timeStep      : no of timestep
; levIndex      : vertical level index
; mapLLC        : (lon,lat) value array of the Lower Left Corner of the map
; mapURC        : (lon,lat) ------- || ------- Upper Right Corner of the map
;                  use -180 to 180 for longitude !
; scaleFactor   : optional scale factor
; selMode       : mode for level spacing of the plot:
;                  halflog (uses half logarythmic levels,i.e. 1,2,5 per decade), 
;                  manual (automatic linear spacing with usage of minVar and maxVar)
;                  auto (default: let do ncl the spacing)
; minVar/maxVar : min/maximal Value to plot for selMode=manual or halflog
; numLevs       : set a number of labels for manual selMode only. NumCol=numLevs+2
;                 (Not set directly, but used for computing the contour level spacing. default=10)
; contLevs      : set individual contour levels (not yet)
; mapType       : projection type (default: lonlat), other: ortho (not compatible with showGrid)
; projSat       : projection type Satellite (default: CylindricalEquidistant, not comp. with "ortho"?)
;                  but compatible with showGrid and mapLLC/URC
; mapLine       : (logical) draws continent lines on the plot (foreground/transparent)
; maskName      : variable to mask with. maskName is expected NOT to have time dimension
;                  'maskName="none"' is same as default (no mask variable)
; maskFile      : optional Filename of the mask variable maskName
;                  it's only taken into account, if maskName is given
; lStrg         : left string
; rStrg         : right string
; tStrg         : title string
; bStrg         : base string - default is prgr name and time stamp only
; showPreview   : (logical) draw a preview window
; maxView       : (logical) maximize plot area on paper (not for buildbot -> convert to png)
; colormap      : string for predefined colormaps of the plot
; showGrid      : display polygon lines with filled colors instead of contour plot (experimental)
; showNcd       : display NDC Grid to find Normalized Device Coordinates on the plot
;---------------------------------------------------------------
; varName=VELOC : test of vector plot - not yet successful - 2011-03
;---------------------------------------------------------------

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;-------------------------------------------------------------------------------
; some global defaults
  DEFAULTLON = "clon"
  DEFAULTLAT = "clat"
  RAD2DEG    = 45./atan(1.)
  PLOTMODE   = "scalar"
  ABORTMSG   = "Abort Script."
  NUMLEVS    = 10
;-------------------------------------------------------------------------------
function sign(val)
begin
  if (val.eq.0) then
    return(1)
  else
    return(val/abs(val))
  end if 
end
;-------------------------------------------------------------------------------
function has_time(dims)
begin
  if (any(dims.eq."time")) then
    return(True)
  else
    return(False)
  end if
end
;-------------------------------------------------------------------------------
function var_has_time(varname,filehandle)
begin
  dims  = getfilevardims(filehandle,varname)
  if ( ismissing(any(dims.eq."time")) ) then
    return(False)
  else
    return(True)
  end if
end
;-------------------------------------------------------------------------------
function selIconVar(varname,filehandle,timestep,levelindex)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 3) then        ; dims: (time,lev,cell)
        var = filehandle->$varname$(timestep,levelindex,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (time,cell)
        var = filehandle->$varname$(timestep,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    else
; TODO: support for 4D tracer fields
;     if (noOfDims .eq. 3) then        ; dims: (lev,cell,tracerID)
;       var = filehandle->$varname$(levIndex,:,tracer)
;       usableVar = True
;     end if
      if (noOfDims .eq. 2) then        ; dims: (lev,cell)
        var = filehandle->$varname$(levelindex,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (cell)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varname + " from input file")
      exit
    end if
    return var
end
;-------------------------------------------------------------------------------
function selRegularVar(varname,filehandle,timestep,levelindex)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
        var = filehandle->$varname$(timestep,levelindex,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 3) then        ; dims: (time,x,y)
        var = filehandle->$varname$(timestep,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    else
      if (noOfDims .eq. 3) then        ; dims: (lev,x,y)
        var = filehandle->$varname$(levelindex,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (x)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varname + " from input file")
      exit
    end if
    return var
end
;-------------------------------------------------------------------------------
procedure scaleVar(var,scalefactor)
local scalefactor
begin
    if (scalefactor .ne. 1) then
      var = var*scalefactor
    end if
end
;-------------------------------------------------------------------------------
function getGridType(variable)
begin
  str   = variable@coordinates
  gtype = stringtocharacter(str)
  if (gtype(0) .eq. "c") then
    str = "cell"
  end if
  if (gtype(0) .eq. "e") then
    str = "edge"
  end if
  if (gtype(0) .eq. "v") then
    str = "vertex"
  end if
  return str
end
;-------------------------------------------------------------------------------
function getCoordinates(variable,filehandle)
begin
  lonlat = str_split(variable@coordinates," ")
  lon    = lonlat(0)
  lat    = lonlat(1)
  x      = filehandle->$lon$
  y      = filehandle->$lat$

  return (/x,y/)
end
;-------------------------------------------------------------------------------
function getLonLats(variable,filehandle)
begin
  lonlat  = str_split(variable@coordinates," ")
  lon     = lonlat(0)
  lat     = lonlat(1)
  x       = filehandle->$lon$*RAD2DEG
  y       = filehandle->$lat$*RAD2DEG

  return (/x,y/)
end
;-------------------------------------------------------------------------------
function createLevels(minVar, maxVar, scaleLimit)
begin
  delimiter = "|"
  tics4positive     = (/1, 2, 5/)
  tics4negative     = (/5, 2, 1/)

  signMin  = sign(minVar)
  signMax  = sign(maxVar)
  if (minVar.eq.0) then 
    logfirst = 0.0
    first = 0.0
  else
    logfirst = ceil(log10(abs(minVar)))
    first = signMin*10^logfirst
  end if
  if (maxVar.eq.0) then
    loglast = 0.0
    last  = 0.0
  else
    loglast  = ceil(log10(abs(maxVar)))
    last  = signMax*10^loglast
  end if


  retval = (/first, last, logfirst, loglast, signMin, signMax, logfirst*signMin - loglast*signMax/)

  myScale    = logfirst
  levels     = ""
  if (first .eq. 0) then
    levels = str_concat((/levels,delimiter,"0"/))
  end if 
  current    = first
  if (current .lt. 0) then
    tics = tics4negative
    do scaleChange=1,scaleLimit
      do ticIndex=0,2
        current = - tics(ticIndex)*10^(logfirst-scaleChange)
        levels = str_concat((/levels,delimiter,flt2string(doubletofloat(current))/))
      end do
    end do
  end if
  if (first.lt.0 .and. last.gt.0) then
    levels = str_concat((/levels,delimiter,"0"/))
  end if
  if (last .gt. 0) then
    tics = tics4positive
    if (first.lt.0) then
      startScale = min((/logfirst, loglast/)) - scaleLimit
    else
      startScale = min((/logfirst, loglast/))
    end if
    j = 0
    do while (current.lt.(last/10))
      myScale = startScale + j
      do ticIndex=0,2
        current = tics(ticIndex)*10^(myScale)
        levels = str_concat((/levels,delimiter,flt2string(doubletofloat(current))/))
      end do
      j = j + 1
    end do
    levels = str_concat((/levels,delimiter,flt2string(doubletofloat(last))/))
  end if

  return(str_split(levels,"|"))
end
;---------------------------------------------------------------
procedure printVarNames(filehandle)
begin
  filevarnames = getfilevarnames(filehandle)
  print("variable names in file: "+filevarnames)
end
;---------------------------------------------------------------
procedure printVar(varname, filehandle)
begin
  dims  = getfilevardims(filehandle,varname)
  sizes = filevardimsizes(filehandle,varname)
  var   = filehandle->$varname$
  ; Print info about variable, its dimensions and attributes
  print(dimsizes(sizes) + " Dimensions:")
  if(.not.any(ismissing(dims))) then
    do j = 0, dimsizes(dims) -1
      print( j + ") " + dims(j) + ": " + sizes(j))
    end do
  end if
  atts = getfilevaratts(filehandle,varname)
  if(.not.any(ismissing(atts))) then
    do k = 0, dimsizes(atts) -1
      print(atts(k) + ": " +filehandle->$varname$@$atts(k)$)
      if (atts(k) .eq. "coordinates") then
        coordinates = var@$atts(k)$
      else
        coordinates = DEFAULTLON + " " + DEFAULTLAT
      end if
    end do
  end if
end
;---------------------------------------------------------------
procedure checkVarname(varname, filehandle)
begin
  filevarnames = getfilevarnames(filehandle)
  if (.not. any(filevarnames.eq.varname)) then
    print("Could not find variable "+varname+"!")
    exit
  end if
  return True
end
;---------------------------------------------------------------
procedure checkDimsOfVars(varname0,varname1,filehandle)
begin
  if (str_join(getfilevardims(filehandle,varname0),"") .ne. str_join(getfilevardims(filehandle,varname1),"")) then
    print("Variables '"+varname0+"' and '"+varname1+"' must have the same dimensions!")
    print(ABORTMSG)
    exit
  end if
end
;---------------------------------------------------------------
function getMaskVar(maskvarname,filehandle,hasOtherMaskfile,timestep,levelindex,plotmode)
begin
  if (hasOtherMaskfile) then
    filehandle = addfile( maskFile, "r" )
  end if
  checkVarname(maskvarname, filehandle)
  if (plotmode.eq."scalar") then
    maskvar = selIconVar(maskvarname,filehandle,timestep,levelindex)
  else
    maskvar = selRegularVar(maskvarname,filehandle,timestep,levelindex)
  end if

  return maskvar
end
;---------------------------------------------------------------
function getMaskDim(maskvarname,filehandle,hasOtherMaskfile)
begin
  if (hasOtherMaskfile) then
    filehandle   = addfile( maskFile, "r" )
  end if
  checkVarname(maskvarname, filehandle)
  maskdims = getfilevardims(filehandle,maskvarname)

  return maskdims
end
;---------------------------------------------------------------
function remapForVecPlot(iFile,resolution,useMask,debug)
begin
  remapFilename="remapnn_"+iFile
  print("#=====================================================================================")
  print("Looking for remapped file: "+remapFilename)
  remapedFiles = systemfunc("ls "+remapFilename)
  if (isstring(remapedFiles(0)) .and. remapedFiles(0).eq.remapFilename) then
    print("Take "+remapFilename+" for accessing vector variables "+str_join(vecVars,","))
    print("Print remove this intermediate file, if it was NOT automatically created from your inpur file "+iFile+"!")
  else
    print("Use CDO to perform remapping: Create intermediate file: "+remapFilename)
    print("Remove this intermediate file, if it was NOT automatically created from your inpur file "+iFile+"!")
    if ( useMask ) then
      cmd = "cdo -P 8 -remapnn,"+resolution+" -selname,"+maskName+","+str_join(vecVars,",")+" "+iFile+" "+remapFilename
      if debug then
        print(cmd)
      end if
      system(cmd)
    else
      cmd = "cdo -P 8 -remapnn,"+resolution+" -selname,"+str_join(vecVars,",")+" "+iFile+" "+remapFilename
      if debug then
        print(cmd)
      end if
    system(cmd)
    end if
  end if
  return remapFilename
end
;---------------------------------------------------------------
procedure checkRemappedFile(infilename,remapfilename)
begin
  orgFile   = addfile (infilename , "r")
  remapFile = addfile(remapfilename,"r")
  orgDims   = getvardims(orgFile)
  remapDims = getvardims(remapFile)
  orgTime   = getfilevardimsizes(orgFile,"time")
  remapTime = getfilevardimsizes(remapFile,"time")
  if ( orgTime .eq. remapTime ) then
    infilename = remapfilename
  else
    print("#=====================================================================================")
    print("Remapped File '"+remapFilename+"' seems to be wrong: Has different number of timesteps")
    print("Original file has "+orgTime(0)+" timesteps")
    print("Remapped File has "+remapTime(0)+" timesteps")
    print(ABORTMSG)
    exit
  end if
end
;---------------------------------------------------------------
function getBoundsOfCoordinates(variable,filehandle)
begin
    lonlat = str_split(variable@coordinates," ")
    lon    = lonlat(0)
    lat    = lonlat(1)


    boundslonName = filehandle->$lon$@bounds
    boundslatName = filehandle->$lat$@bounds
    boundslon     = filehandle->$boundslonName$ * RAD2DEG
    boundslat     = filehandle->$boundslatName$ * RAD2DEG

    return (/boundslon, boundslat/)
end
;---------------------------------------------------------------
procedure setBoundsEnds(blon,blat,lonmin,lonmax,latmin,latmax)
begin
  latmax=max(blat)
  latmin=min(blat)
  lonmax=max(blon)
  lonmin=min(blon)

  if ( lonmin.lt.0 ) then
  ; longitudes are given in the range [-180,180]
    lonmin = -180.
    lonmax =  180.
  else
  ; longitudes are given in the range [0,360]
    lonmin =    0.
    lonmax =  360.
  end if
end
;---------------------------------------------------------------
procedure checkLongitude(x,boundslon,lonmin,lonmax,debug)
begin
  dlon      = 80.
  ncell_tot = dimsizes(x)  ; total # of cells
  nfix      = 0            ; count the # of problematic cells

  do icell = 0, ncell_tot-1
    if ( any(boundslon(icell,:).le.(lonmin+dlon)) .and. \
         any(boundslon(icell,:).gt.(lonmax-dlon))       ) then

      bl1=boundslon(icell,0)
      bl2=boundslon(icell,1)
      bl3=boundslon(icell,2)

      boundslon(icell,:) = where(boundslon(icell,:).gt.(lonmax-dlon), \;
                           boundslon(icell,:)-360.,             \; where true
                           boundslon(icell,:) )                  ; where false

      bn1=boundslon(icell,0)
      bn2=boundslon(icell,1)
      bn3=boundslon(icell,2)

      if (debug) then
        print("cell="+icell+" lon="+bl1+", "+bl2+", "+bl3+" changed into " \
                           +" lon="+bn1+", "+bn2+", "+bn3)
      end if

      nfix = nfix +1
    end if
  end do
  print("Longitude(s) of "+nfix+" cells corrected to avoid plotting problem.")
end
;---------------------------------------------------------------
procedure setLevels(selmode,resource,minvar,maxvar,scalelimit,numlevs,debug)
begin
  if (selmode .eq. "manual") then
    resource@cnLevelSelectionMode = "ManualLevels"
    resource@cnMinLevelValF       = minvar
    resource@cnMaxLevelValF       = maxvar
;    diffLog10                 = log10(abs(maxVar-minVar))
;;   if (diffLog10 .lt. 0) then
;      resource@cnLevelSpacingF    = 10^(floor(diffLog10))/numLevs
;;   else
;;     resource@cnLevelSpacingF    = 10^(floor(diffLog10))/numLevs
;;   end if
    diffspacing               = abs(maxvar-minvar)/(numlevs)
    resource@cnLevelSpacingF      = diffspacing
    if (debug) then
      print("Manual Plotlevel setting......")
      print("  minVal  : " + minvar)
      print("  maxVal  : " + maxvar)
      print("  numLevs: " + numlevs)
    end if
  end if

  if (selmode .eq. "halflog") then
    resource@cnLevelSelectionMode = "ExplicitLevels"
    plotLevels = createLevels(minvar,maxvar,scalelimit)
    if (debug) print("plotlevels = "+plotLevels) end if
    resource@cnLevels             = plotLevels
  end if
end
;---------------------------------------------------------------
procedure setDefaultResource(resource)
begin
  resource@gsnMaximize          = False
; resource@gsnPaperOrientation  = "landscape"  ; turns plot on paper
; resource@wkPaperSize          = "A4"
  resource@gsnFrame             = False
  resource@gsnSpreadColors      = True

  resource@cnFillOn             = True   ;  writes labelbar, no contour labels
  resource@cnLinesOn            = True
  resource@cnInfoLabelOn        = False

  FontHeight0               = 0.012
  FontHeight1               = 0.015
  FontHeight2               = 0.017
  FontHeight3               = 0.020

  resource@tiXAxisFontHeightF   = FontHeight0
  resource@tiYAxisFontHeightF   = FontHeight0
  resource@tmXBLabelFontHeightF = FontHeight0
  resource@tmYLLabelFontHeightF = FontHeight0
  resource@tmXBLabelJust        = "CenterCenter"
  resource@gsnStringFontHeightF = FontHeight0

  resource@mpFillOn               = True
  resource@cnFillOn               = True

  resource@lbLabelBarOn             = True
  resource@pmLabelBarHeightF        = 0.12
  resource@pmLabelBarWidthF         = 0.80
  resource@pmLabelBarOrthogonalPosF = 0.15
  resource@lbLabelAutoStride        = True
  resource@lbLabelFontHeightF       = FontHeight1     ;   color bar labels

  resource@cnFillMode           = "rasterfill"
  resource@cnRasterSmoothingOn  = True
  resource@mpGreatCircleLinesOn = True

end
;---------------------------------------------------------------
procedure setPlotCaptions(resource,leftstring,rightstring,centerstring,titlestring)
begin
  resource@gsnLeftString   = leftstring
  resource@gsnRightString  = rightstring
  resource@gsnCenterString = centerstring
  vesource@tiMainString    = titlestring
end
;---------------------------------------------------------------
procedure setAutomaticPlotCaptions(resource,plotmode,varname,filehandle,filename,timestep,levelindex)
local varname
begin
  ; titles
  if(.not. isvar("lStrg")) then
    hastime = var_has_time(varname, filehandle)
    if (hastime) then
      timeval = filehandle->time(timestep)
      dayfrag = timeval - floor(timeval)
      day     = doubletointeger(floor(timeval))
      hour    = doubletointeger(floor(86400*dayfrag/3600))
      minute  = doubletointeger(floor((86400*dayfrag - hour*3600)/60))
      second  = doubletointeger(86400*dayfrag - hour*3600 - minute*60)
      timeStr = str_join((/day,hour,minute,second/),"|")
      if (plotmode.eq."vector") then
        varname = "Velocity"
      end if
      lStrg   = (/varname,"(plstep "+timestep+"["+timeStr+"],","pllev "+levelindex+")"/)
   else
     lStrg   = (/varname,"(","pllev "+levelindex+")"/)
   end if
  end if
  leftstring   = str_join(lStrg," ")

  if(.not. isvar("rStrg")) then
    rStrg      = (/filename/)
  end if
  rightstring  = str_join(rStrg," ")

  if(.not. isvar("tStrg")) then
    tStrg        = "ICON"
  end if
  titlestring  = str_join(tStrg," ")
  centerstring = ""

  resource@gsnLeftString   = leftstring
  resource@gsnRightString  = rightstring
  resource@gsnCenterString = centerstring
  resource@tiMainString    = titlestring
end
;---------------------------------------------------------------
procedure setMapType(resource,maptype)
begin
  if (maptype .eq. "ortho") then
    resource@mpProjection          = "Orthographic"
    resource@mpPerimOn             = False    ; turn off box around plot
    resource@mpGridAndLimbOn       = True     ; draw grid lines and limb line
    resource@mpGridLineDashPattern = 2        ; choose pattern used to draw the grid

  else
    resource@mpProjection         = "CylindricalEquidistant"
    resource@mpLimitMode          = "LatLon"
  end if
end
;---------------------------------------------------------------
procedure selMapCut(resource,mapllc,mapurc)
begin
  ; bounds of plotting area (not used for Orthographic)
  resource@mpMinLatF    = mapllc(1)
  resource@mpMaxLatF    = mapurc(1)
  resource@mpMinLonF    = mapllc(0)
  resource@mpMaxLonF    = mapurc(0)

  ; center of view (for Orthographic and Sattelite only)
  resource@mpCenterLonF = (resource@mpMaxLonF+resource@mpMinLonF)/2
  resource@mpCenterLatF = (resource@mpMaxLatF+resource@mpMinLatF)/2
end
;---------------------------------------------------------------
procedure  setSateliteProj(resource,projsat)
begin
  if (projsat) then
    resource@mpProjection = "Satellite"
    resource@mpMinLatF    = resource@mpMinLatF + 10
    resource@mpMaxLatF    = resource@mpMaxLatF - 10
    resource@mpMinLonF    = resource@mpMinLonF + 10
    resource@mpMaxLonF    = resource@mpMaxLonF - 10
  end if
end
;---------------------------------------------------------------
procedure setMapVisibility(resource,mapline)
begin
  if (mapline) then
    resource@mpGeophysicalLineColor = "foreground"
  else
    resource@mpGeophysicalLineColor = "transparent"
  end if
end
;---------------------------------------------------------------
procedure checkMinMaxVar(minvar,maxvar)
begin
  if ( minvar .gt. maxvar ) then
    print("minVar has to be larger than maxVar")
    exit
  end if
end
;---------------------------------------------------------------
function getLshiftForNDCString(mystring)
begin
  lshift = (strlen(mystring)/2 - 1)/100.0
  if (strlen(mystring) .gt. 4)
  lshift = lshift - lshift/5 + 0.1/strlen(mystring)
  else
  lshift = lshift + 0.01
  end if

  return lshift
end
;---------------------------------------------------------------
procedure setBaseString(workstation,plot,resource,basestring)
begin
  resource@txFontHeightF = 0.01 ; if this is changed, please adjust the implementation of getLshiftForNDCString()
  gsn_text_ndc(workstation,basestring,getLshiftForNDCString(basestring),.01,resource)
end
;---------------------------------------------------------------
procedure setAutomaticBaseString(workstation,resource)
begin
  if(.not. isvar("bStrg")) then
    bStrg      = ""
  end if
  basestr     = str_join(bStrg," ")
  ; base string: plot file name + basestring + time stamp
  ; attention: virtual page changes when lon/lat ratio changes
  BaseString = "Prgr icon_ocean.ncl: " + " " + basestr + " " + systemfunc("date")
  resource@txFontHeightF = 0.011
  gsn_text_ndc(workstation,BaseString,getLshiftForNDCString(BaseString),.12,resource)  ; TODO: string align left
end
;---------------------------------------------------------------
procedure setMaskColor(wks,resource)
begin
  newcolor = NhlNewColor(wks,0.85,0.85,0.85)   ; add gray to colormap
  ;resource@gsnSpreadColorStart  = 1
  resource@gsnSpreadColorEnd     = 19  ; last color number not known?
  resource@cnMissingValFillColor = 20  ; last color is used for missing values
end
;---------------------------------------------------------------
procedure setCoordinates(resource,lons,lats)
begin
  resource@sfXArray = lons
  resource@sfYArray = lats
end
;---------------------------------------------------------------
procedure showMapInfo(resource,maptype,projsat,mapline)
begin
  print("mapType: "+maptype)
  print("projSat: "+projsat)
  print("mapLine: "+mapline)
  print("CentLon/Lat="+resource@mpCenterLonF+"; "+resource@mpCenterLatF)
  print("Min/MaxLatF="+resource@mpMinLatF+"; "+resource@mpMaxLatF)
  print("Min/MaxLonF="+resource@mpMinLonF+"; "+resource@mpMaxLonF)
end
;---------------------------------------------------------------
function getFlags(var,boundslon,boundslat,resource)
begin
  flags = new(dimsizes(var),logical,"No_FillValue")
  do i = 0,dimsizes(var) - 1
    flags(i) = where(all(boundslon(i,:) .gt. resource@mpMaxLonF) .or. \
                     all(boundslon(i,:) .lt. resource@mpMinLonF) .or. \
                     all(boundslat(i,:) .gt. resource@mpMaxLatF) .or. \
                     all(boundslat(i,:) .lt. resource@mpMinLatF), \
                     False, True)
  end do
  return flags
end
;---------------------------------------------------------------
procedure setDefaultVectorPlot(resource, refmagnitude, reflength, style, mindistance)
begin
  resource@vcRefMagnitudeF         = refmagnitude ; make vectors larger
  resource@vcRefLengthF            = reflength    ; reference vector length
  resource@vcGlyphStyle            = style        ; turn on curly vectors
  resource@vcMinDistanceF          = mindistance  ; thin out vectors
end
;---------------------------------------------------------------
function getColorsFromPlot(plot)
begin
  getvalues plot@contour
    "cnFillColors" : colors
  end getvalues
  return colors
end
;---------------------------------------------------------------
function getLevelsFromPlot(plot)
begin
  getvalues plot@contour
    "cnLevels"     : levels
  end getvalues
  return levels
end
;---------------------------------------------------------------
procedure printGridPlotInfo(colors,levels,flags)
begin
  print(colors)
  print(levels)
  print ("Outside the plot area:           " + dimsizes(ind(flags .eq. False)) + \
      " triangles - not plotted")
end
;---------------------------------------------------------------
function setupGridResource(var,levels,colors,boundslon,boundslat,wks,plot,debug)
begin
  pres             = True
  pres@gsEdgesOn   = True   ; Turn on edges
  pres@gsFillIndex = 0      ; Solid fill, the default

  ; First draw the triangles associated with the lowest level.
  i = 0
  vlow = ind(var .lt. levels(i))
  ; if no index vlow is found with values less than levels(0), vlow is missing value:
  if (.not. ismissing(vlow(0))) then
    do j = 0, dimsizes(vlow)-1
      pres@gsFillColor = colors(i)   ; first color
      gsn_polygon(wks,plot,boundslon(vlow(j),:),boundslat(vlow(j),:),pres)
    end do
  end if
  if (debug) then
    print("i=" + i + ", Values smaller than       "+levels(i)+": " + dimsizes(vlow) + \
          " triangles considered, color= " + colors(i))
  end if
  delete(vlow)

  ; Now draw the triangles between the lowest and highest levels.
  do i = 0, dimsizes(levels) -2
    vind = ind(var .ge. levels(i) .and. var .lt. levels(i+1))
    ; if no index is found, vind is missing value:
    if (.not. ismissing(vind(0))) then
    do j = 0, dimsizes(vind)-1
      pres@gsFillColor = colors(i+1)
                       gsn_polygon( wks,plot,   boundslon(vind(j),:),boundslat(vind(j),:),pres)
    end do
    end if
    if (debug) then
      print ("i=" + i + ", Values between " + levels(i) + " and " + levels(i+1) +": "\
             +dimsizes(vind) + " triangles considered, color= " + colors(i+1))
    end if
    delete(vind)
  end do

  ; Finally draw the triangles associated with the highest level.
  i = dimsizes(levels) -1
  vhig = ind(var .ge. levels(i) )

  ; if no index vhig is found with values larger levels(max-1), vhig is missing value:
  if (.not. ismissing(vhig(0))) then
  do j = 0, dimsizes(vhig) -1
    pres@gsFillColor = colors(i+1)
                     gsn_polygon( wks,plot,   boundslon(vhig(j),:),boundslat(vhig(j),:),pres)
  end do
  end if
  if (debug) then
    print ("i=" + i + ", Values   larger than      "+levels(i) +": " \
           + dimsizes(vhig) + " triangles considered, color= " + colors(i+1))
  end if
  delete(vhig)

  return pres
end
;---------------------------------------------------------------
procedure plotGrid(workstation,plot,variable,x,filehandle,resource,debug)
begin
  bounds    = getBoundsOfCoordinates(variable,filehandle)
  boundslon = bounds(0,:,:)
  boundslat = bounds(1,:,:)

  latmax=flt2dble(0.0)
  latmin=flt2dble(0.0)
  lonmax=flt2dble(0.0)
  lonmin=flt2dble(0.0)
  setBoundsEnds(boundslon,boundslat,lonmin,lonmax,latmin,latmax)
  checkLongitude(x,boundslon,lonmin,lonmax,debug)

  flags = getFlags(variable,boundslon,boundslat,resource)

  colors = getColorsFromPlot(plot)
  levels = getLevelsFromPlot(plot)

  if (debug) printGridPlotInfo(colors,levels,flags) end if

  pres = setupGridResource(variable,levels,colors,boundslon,boundslat,workstation,plot,debug)
end
;===============================================================
;----- MAIN PROGRAM --------------------------------------------
 begin
;- options handling + checking
  plotMode  = PLOTMODE

  wcStrt    = systemfunc("date")

  print("+++++ "+wcStrt)

  Model     = "Icon: explicit ocean test"

  ; input/output files
  if(.not. isvar("iFile")) then
    print("Input file is required! Use iFile option.")
    print(ABORTMSG)
    exit
  else
    if ( .not. isfilepresent(iFile) )
      print("Could not read from input file: "+iFile+"!")
      print(ABORTMSG)
      exit
    end if
  end if
  if(.not. isvar("oType")) oType = "eps" end if
  if(.not. isvar("oFile")) then
    ext      = get_file_suffix(iFile,0)
    oFile    = ext@fBase
  end if
  print("Outputfile '"+str_concat((/oFile,".",oType/))+"' will be creates in "+systemfunc("dirname "+iFile))

  if(.not. isvar("numLevs")) then
    numLevs  = NUMLEVS
  else
    if (numLevs .lt. 1) then
      print("WARNING: numLevs must be >= 1")
      print(ABORTMSG)
      exit
    end if
  end if

  if (.not. isvar("DEBUG"))       DEBUG       = False           end if
  if (.not. isvar("showPreview")) showPreview = False           end if
  if (.not. isvar("showGrid"))    showGrid    = False           end if
  if (.not. isvar("showNdc"))     showNdc     = False           end if
  if (.not. isvar("maxView"))     maxView     = False           end if
  if (.not. isvar("scaleLimit"))  scaleLimit  = 2               end if
  if (.not. isvar("selMode"))     selMode     = "auto"          end if
  if (.not. isvar("varName"))     varName     = "ELEV"          end if
  if (.not. isvar("timeStep"))    timeStep    = 0               end if
  if (.not. isvar("levIndex"))    levIndex    = 0               end if
  if (.not. isvar("scaleFactor")) scaleFactor = 1               end if
  if (.not. isvar("colormap"))    colormap    = "BlueDarkRed18" end if

  if (.not.isvar("mapLLC"))       mapLLC      = (/-180.,-90./)  end if
  if (.not.isvar("mapURC"))       mapURC      = (/180.0,90.0/)  end if
  if (.not.isvar("mapType"))      mapType     = "lonlat"        end if
  if (.not.isvar("projSat"))      projSat     = False           end if
  if (.not.isvar("mapLine"))      mapLine     = True            end if
  if (.not.isvar("resolution"))   resolution  = "r90x45"        end if
  ; default + optional coordinates
  if (isvar("lonCo"))             DEFAULTLON  = lonCo           end if
  if (isvar("latCo"))             DEFAULTLAT  = latCo           end if

  ; mask handling
  if (.not.isvar("maskName"))
    useMask     = False
  else
    if ( .not.isstring(maskName)) then
      print("Parameter 'maskName' has to be a string!")
      exit
    else
      useMask  = True
      if (.not.isvar("maskFile"))
        if (maskName.eq."none") then
          useMask = False
        else
          if (DEBUG) then
            print("Use internal mask variable " +maskName)
          end if
        end if
      else
        if (DEBUG) then
          print("Use external mask file " + maskFile)
        end if
      end if
    end if
  end if


  ; expect a regular grid for vector plot as you get it from
  if (isvar("vecVars")) then
    ; expectes is a string or an array of size 2: "u-veloc v-veloc" or (/"u-veloc","v-veloc"/)
    if (isstring(vecVars)) then
      vecVars_ = str_split(vecVars," ")
      delete(vecVars)
      vecVars = vecVars_
    end if
    plotMode = "vector"
  end if

  if (plotMode.eq."vector") then
    ; performe some remapping to a regular grid because ncl cannot draw vector
    ; from unstructured grids
    remapFilename = remapForVecPlot(iFile,resolution,useMask,DEBUG)
    checkRemappedFile(iFile,remapFilename)
  end if

  File = addfile( iFile, "r" )

  if (DEBUG) printVarNames(File) end if

  ;Read mask variable
  if( useMask ) then
    maskVar = getMaskVar(maskName,File,isvar("maskFile"),timeStep,levIndex,plotMode)
    maskDim = getMaskDim(maskName,File,isvar("maskFile"))
  end if

  if (plotMode.eq."vector") then
    varnames = vecVars
  else
    varnames = (/varName/)
  end if
  print("Plot variables: " + varnames)

  ; Reading data variables =====================================================
  if (plotMode.eq."scalar") then ; scalar mode ================================

    printVar(varName, File)

    var = selIconVar(varName,File,timeStep,levIndex)

    scaleVar(var,scaleFactor)

    lonlats = getLonLats(var,File)
    x       = lonlats(0,:)
    y       = lonlats(1,:)

    if (DEBUG) if (useMask) print(maskDim) end if end if

    ; set minVar_maxVar for plotting
    if(.not. isvar("minVar")) minVar = min(var) end if
    if(.not. isvar("maxVar")) maxVar = max(var) end if
    checkMinMaxVar(minVar,maxVar)

    if ( useMask ) then
      ; set variable var to missing, where var is not equal mvalue (3rd
      ; parameter)
      var = mask(var,maskVar,1)
      ;slm = maskVar - 0.5
      ;var = mask ( var, slm, 0.5)
    end if

  else ; vector mode ==========================================================

    uvarname = varnames(0)
    vvarname = varnames(1)

    checkDimsOfVars(uvarname,vvarname,File)

    uvar = selRegularVar(uvarname,File,timeStep,levIndex)
    vvar = selRegularVar(vvarname,File,timeStep,levIndex)

    scaleVar(uvar,scaleFactor)
    scaleVar(vvar,scaleFactor)

    if ( useMask ) then
      uvar = mask(uvar,maskVar,1)
      vvar = mask(vvar,maskVar,1)
    end if

    velocity = sqrt(uvar*uvar + vvar*vvar)
    ; set minVar/maxVar for plotting
    if(.not. isvar("minVar")) minVar = min(velocity) end if
    if(.not. isvar("maxVar")) maxVar = max(velocity) end if
    checkMinMaxVar(minVar,maxVar)
  end if ; Reading data variables =============================================

  ;---------------------------------------------------------------
  if (DEBUG) then
    if ( var_has_time(varName, File) ) then
      timeVal = File->time(timeStep)
      print("time "+timeVal)
    end if
    print("iFile       = "+iFile)
    print("oFile       = "+oFile)
    print("Graphics format is " +oType)
    print("varName     = "+varName)
    print("timeStep    = "+timeStep)
    print("mapLLC(lon) = "+mapLLC(0))
    print("mapLLC(lat) = "+mapLLC(1))
    print("mapURC(lon) = "+mapURC(0))
    print("mapURC(lat) = "+mapURC(1))
    print("minVar      = "+minVar)
    print("maxVar      = "+maxVar)
    if (useMask) then
    print("maskName    = "+maskName)
    end if
  end if
  ;---------------------------------------------------------------

  ;---------------------------------------------------------------
  ; make the plot
  ;---------------------------------------------------------------
  ; preparations
  wks = gsn_open_wks(oType,oFile)

  if (showNdc) drawNDCGrid(wks) end if

  gsn_define_colormap(wks,colormap)

  if (showPreview) then
    xwks = gsn_open_wks("x11","x11")
    gsn_define_colormap(xwks,colormap)
  end if

  ResC = True

  setDefaultResource(ResC)

  if (maxView) ResC@gsnMaximize = True end if

  if (useMask) setMaskColor(wks,ResC) end if

  if (plotMode.eq."scalar") setCoordinates(ResC,x,y) end if

  setAutomaticPlotCaptions(ResC,plotMode,varName,File,iFile,timeStep,levIndex)

  setAutomaticBaseString(wks,ResC)
  if (showPreview) setAutomaticBaseString(xwks,ResC) end if

  setLevels(selMode,ResC,minVar,maxVar,scaleLimit,numLevs,DEBUG)

  setMapType(ResC,mapType)

  selMapCut(ResC,mapLLC,mapURC)

  setSateliteProj(ResC,projSat)

  setMapVisibility(ResC,mapLine)

  if (DEBUG) print("Gridtype is "+getGridType(var)) end if
  if (DEBUG) showMapInfo(ResC,mapType,projSat,mapLine) end if

  ; MAIN PLOT CALLS ===========================================================
  if (plotMode .eq. "vector") then ; vector plot ==============================
    if (showGrid) then
      print("#= WARNING =============================================")
      print("Display Vectors and the underlying grid is not usefull, ")
      print("because original data is interpolated to a regular grid for vector representation.")
      print(ABORTMSG)
      exit
    else
      setDefaultVectorPlot(ResC,5.0,8.0,"CurlyVector",0.017)

      if ( showPreview ) then
        x11 = gsn_csm_vector_scalar_map(xwks,uvar,vvar,uvar*uvar+vvar*vvar,ResC)
      end if
      vc = gsn_csm_vector_scalar_map(wks,uvar,vvar,uvar*uvar+vvar*vvar,ResC)
    end if
  else ; contour plot =========================================================
    plot                     = gsn_csm_contour_map(wks,var,ResC)
    if (showPreview) preview = gsn_csm_contour_map(xwks,var,ResC) end if

    if (showGrid) then
      if (showPreview) then
        plotGrid(xwks,preview,var,x,File,ResC,DEBUG)
        print("#=======================================================")
        print("For reasons of performance, the preview of a grid plot is not written to the output file be default")
        print("Please set 'showPreview=False' on the command line!")
      else
        plotGrid(wks,plot,var,x,File,ResC,DEBUG)
      end if
    end if
  end if

  if (showPreview) frame(xwks) end if
  frame(wks)
end
