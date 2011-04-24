;---------------------------------------------------------------
; Usage:
;   ncl icon_ocean.ncl 'iFile="oceanTest/testfile.nc"' 'varName="W"' minVar=-1.0 maxVar=0.5 'mapLLC=(/35.0, -8/)'  'mapURC=(/55, 8/)' timeStep=1 levIndex=1
;   ncl icon_ocean.ncl 'iFile="oceanTest/testfile.nc"' 'varName="ELEV"' minVar=-1.0 maxVar=0.5 'mapLLC=(/35.0, -8/)'  'mapURC=(/55, 8/)' timeStep=1

;   masking:
;     ncl icon_ocean.ncl 'iFile="oce_aqua_planet_O.R2B04L4_0001.nc"' 'varName="ELEV"' 'oFile="test"' 'maskName="wet_c"' 
;     ncl icon_ocean.ncl 'iFile="oce_aqua_planet_O.R2B04L4_0001.nc"' 'varName="ELEV"' 'oFile="test"' 'maskName="topo"' 'maskFile="icon_mask.nc"' 
;
;  vectorplot:
;     ncl icon_ocean.ncl 'iFile="oce.nc"' 'vecVars="u-veloc v-veloc"' 'oFile="test"'
;---------------------------------------------------------------
; This script makes contour plots of ocean state
; Parameter:
; iFile         : input file
; oFile         : plot file
; oType         : output graphics format (ps, eps, png, ...)
; varName       : name of the variable to (contour) plot
; vecVars       : space separated string of the 2 vector components u and v to draw a vector plot
; resolution    : resolution string whish is used for remapping the icon data to regular grid
; minVar/maxVar : min/maximal Value to plot
; timeStep      : no of timestep
; levIndex      : level index
; mapLLC        : (lon,lat) value array of the Lower Left Corner of the map
; mapURC        : (lon,lat) ------- || ------- Upper Right Corner of the map
;                  use -180 to 180 for longitude !
; scaleFactor   : optional scale factor
; selMode       : mode for level spacing of the plot:
;                  halflog (uses half logarythmic levels,i.e. 1,2,5 per decade), 
;                  manual (automatic linear spacing with usage of minVar and maxVar)
;                  auto (default: let do ncl the spacing)
; contLevs      : set a number of labels for manual selMode only. It will not set this number
;                 directly, but is used for computing the contour level spacing. default=10
; mapType       : projection type (default: lonlat), other: ortho (not compatible with showGrid)
; maskName      : variable to mask with. maskName is expected NOT to have time dimension
; maskFile      : optional Filename of the mask variable maskName - it's only taken into account, if maskName is given
; lStrg         : left string
; rstrg         : right string
; tStrg         : title string
; bStrg         : base string - default is the time stamp
; showPreview   : (logical) draw a preview window
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
function sign(val)
begin
  if (val.eq.0) then
    return(1)
  else
    return(val/abs(val))
  end if 
end
;-------------------------------------------------------------------------------
function createLevels(minVar, maxVar, minScale)
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
  scaleLimit = 2
  if (first .eq. 0) then
    levels = str_concat((/levels,delimiter,"0"/))
  end if 
  current    = first
  if (current .lt. 0) then
    tics = tics4negative
    do scaleChange=1,scaleLimit
      do ticIndex=0,2
        current = - tics(ticIndex)*10^(logfirst-scaleChange)
        levels = str_concat((/levels,delimiter,flt2string(current)/))
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
       levels = str_concat((/levels,delimiter,flt2string(current)/))
     end do
     j = j + 1
   end do
   levels = str_concat((/levels,delimiter,flt2string(last)/))
 end if

  return(str_split(levels,"|"))
end
;procedure initOptions(plotMode,
;---------------------------------------------------------------
 begin

  plotMode  = "scalar"

  rad2deg   = 45./atan(1.)

  wcStrt    = systemfunc("date")

  print("+++++ "+wcStrt)

  Model     = "Icon: explicit ocean test"

  if(.not. isvar("DEBUG")) then
    DEBUG       = False
  end if
  if(.not. isvar("showPreview")) then
    showPreview = False
  end if
  if(.not. isvar("iFile")) then
    iFile       = "test_hom_1_O.R2B04L4_0001.nc"
  end if
  if(.not. isvar("oFile")) then
    ext      = get_file_suffix(iFile,0)
    oFile    = ext@fBase
   ; ext      = get_file_suffix(systemfunc("basename "+iFile),0)
   ; oFile    = ext@fBase
  end if
  if(.not. isvar("oType")) then
    oType    = "eps"
  end if
  if(.not. isvar("contLevs")) then
    contLevs  = 10
  else
    if (contLevs .lt. 1) then
      print("WARNING: contLevs must be >= 1")
      exit
    end if
  end if
  if(.not. isvar("minScale")) then
    minScale = -1
  end if

  if(.not. isvar("selMode")) then
    selMode     = "auto"
  end if

  if (.not. isvar("varName")) then
    varName     = "ELEV"
  end if
  if (.not. isvar("timeStep")) then
    timeStep    = 0
  end if
  if (.not. isvar("levIndex")) then
    levIndex    = 0
  end if

  if (.not. isvar("scaleFactor")) then
    scaleFactor       = 1
  end if

  if(.not. isvar("colormap")) then
    colormap = "BlueDarkRed18"
  end if

; map info
  if (.not.isvar("mapLLC"))
    mapLLC      = (/-180.,-90./)
  end if
  if (.not.isvar("mapURC"))
    mapURC      = (/180.0,90.0/)
  end if
  if (.not.isvar("mapType"))
    mapType     = "lonlat"
  end if
  if (.not.isvar("maskName"))
    useMask     = False
  else
    if ( .not.isstring(maskName)) then
      print("Parameter 'maskName' has to be a string!")
      exit
    else
      useMask  = True
      if (.not.isvar("maskFile"))
        if (DEBUG) then
          print("Use internal mask variable " +maskName)
        end if
      else
        if (DEBUG) then
          print("Use external mask file " + maskFile)
        end if
      end if
    end if
  end if
  if (.not.isvar("showGrid"))
    showGrid    = False
  end if
  if (.not.isvar("showNdc"))
    showNdc     = False
  end if

  ; default + optional coordinates
  defaultLon = "clon"
  defaultLat = "clat"
  if (isvar("lonCo"))
    defaultLon    = lonCo
  end if
  if (isvar("latCo"))
    defaultLat    = latCo
  end if

  if (.not.isvar("resolution")) then
    resolution="r90x45"
  end if

  ; expect a regular grid for vector plot as you get it from
  ;   cdo -P 8 remapnn,r360x180 -selname,wet_c,u-veloc,v-veloc icon_oce.nc vectors.nc
  if (isvar("vecVars")) then
    ; expectes is a string or an array of size 2: "u-veloc v-veloc" or (/"u-veloc","v-veloc"/)
    if (isstring(vecVars)) then
      vecVars_ = str_split(vecVars," ")
      delete(vecVars)
      vecVars = vecVars_
    end if
    plotMode = "vector"
    ; performe the remapping
    remapFilename="remapnn_"+iFile
    remapedFiles = systemfunc("ls "+remapFilename)
    if (isstring(remapedFiles(0)) .and. remapedFiles(0).eq.remapFilename) then
      print("Take "+remapFilename+" for accessing vector variables "+str_join(vecVars,","))
    else
      print("Use CDO to perform remapping")
      if ( useMask ) then
        cmd = "cdo -P 8 -remapnn,"+resolution+" -selname,"+maskName+","+str_join(vecVars,",")+" "+iFile+" "+remapFilename
        if DEBUG then
          print(cmd)
        end if
        system(cmd)
      else
        cmd = "cdo -P 8 -remapnn,"+resolution+" -selname,"+str_join(vecVars,",")+" "+iFile+" "+remapFilename
        if DEBUG then
          print(cmd)
        end if
      system(cmd)
      end if
    end if
  end if

  if (plotMode.eq."vector") then
   iFile = remapFilename
  end if

  File = addfile( iFile, "r" )

  filevarnames = getfilevarnames(File)
  ; check if varName can be found in iFile
  if (plotMode.eq."vector") then
    do ivar = 0, 1
      if ( .not. any(filevarnames.eq.vecVars(ivar)) ) then
        print("NOPE! Vector variable "+vecVars(ivar)+" could not be found.")
        exit
      end if
    end do
  else
    if ( .not. any(filevarnames.eq.varName) ) then
      print("NOPE! Variable could not be found.")
      exit
    end if
  end if
  ; Check for presence of mask variable
  if( useMask ) then
    if ( .not.isvar("maskFile") ) then
      if ( .not. any(filevarnames.eq.maskName) ) then
        print("NOPE! Mask-Variable could not be found.")
        exit
      end if
    else
      mFile = addfile( maskFile, "r" )
      maskfilevarnames = getfilevarnames(mFile)
      if ( .not. any(maskfilevarnames.eq.maskName) ) then
        print("NOPE! Mask-Variable could not be found.")
        exit
      end if
    end if
  end if

  if (any(filevarnames.eq."time")) then
    timeVal = File->time(timeStep)
  else
    print("Variable has only one dimension! Time is not avaliable")
  end if

  if (plotMode.eq."vector") then 
    varnames = vecVars
  else
    varnames = (/varName/)
  end if
  print("Plot variables: " + varnames)

  ;varArray = new(1,typeof(File->$varnames(0)$),"No_FillValue")

  if (plotMode.eq."scalar") then
    dims  = getfilevardims(File,varName)
    sizes = filevardimsizes(File,varName)
    var   = File->$varName$
    print(dimsizes(sizes) + " Dimensions:")
    if(.not.any(ismissing(dims))) then
      do j = 0, dimsizes(dims) -1
        print( j + ") " + dims(j) + ": " + sizes(j))
      end do
    end if
    atts = getfilevaratts(File,varName)
    if(.not.any(ismissing(atts))) then
      do k = 0, dimsizes(atts) -1
        print(atts(k) + ": " +File->$varName$@$atts(k)$)
        if (atts(k) .eq. "coordinates") then
          coordinates = var@$atts(k)$
        else
          coordinates = defaultLon + " " + defaultLat
        end if
      end do
    end if
    delete(var)

    ;Read mask variable
    if( useMask ) then
      if (isvar("maskFile")) then
        maskVar = mFile->$maskName$
        mdims  = getfilevardims(mFile,maskName)
      else
        maskVar = File->$maskName$
        mdims  = getfilevardims(File,maskName)
      end if
      print("Found maskVar: "+maskName)
      ; Select the correct level
      maskDimsSize = dimsizes(mdims)
      print(mdims)
      print(dims)
      if (2.eq.maskDimsSize) then
        maskVar_ = maskVar(levIndex,:)
        delete(maskVar)
        maskVar = maskVar_
      end if
    end if

  ; Read variable and its coordinates from file
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if (noOfDims .eq. 3) then        ; dims: (time,lev,cell) 
      var = File->$varName$(timeStep,levIndex,:)
      usableVar = True
      checkgrid = dims(2)
    end if
    if (noOfDims .eq. 2) then        ; dims: (time,cell)
      var = File->$varName$(timeStep,:)
      usableVar = True
      checkgrid = dims(1)
    end if
    if (noOfDims .eq. 1) then        ; dims: (cell)
      var = File->$varName$(:)
      usableVar = True
      checkgrid = dims(0)
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varName + " from input file")
      exit
    end if

    if (scaleFactor .ne. 1) then
      var = var*scaleFactor
    end if

    lonlat = str_split(coordinates," ")
    lon    = lonlat(0)
    lat    = lonlat(1)
    x      = File->$lon$ *rad2deg
    y      = File->$lat$ *rad2deg

    if (DEBUG) then
      print(dims)
      print(atts)
      print(coordinates)
    end if

    ; set minVar/maxVar for plotting
    if(.not. isvar("minVar")) then
      minVar = min(var)
    end if
    if(.not. isvar("maxVar")) then
      maxVar = max(var)
    end if
    if ( minVar .gt. maxVar ) then
      print("minVar has to be larger than maxVar")
      exit
    end if

    if ( useMask ) then
      var = mask(var,maskVar,0)
    end if
  else
    uvarname = varnames(0)
    vvarname = varnames(1)
    udims  = getfilevardims(File,uvarname)
    vdims  = getfilevardims(File,vvarname)
    if (str_join(udims,"").ne.str_join(vdims,"")) then
      print("u-Variable and v-Variable must have the same dimensions!")
      exit
    end if
    noOfDims  = dimsizes(udims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
      uvar = File->$uvarname$(timeStep,levIndex,:,:)
      vvar = File->$vvarname$(timeStep,levIndex,:,:)
    end if
    if (noOfDims .eq. 3) then        ; dims: (time,x,y)
      uvar = File->$uvarname$(timeStep,:,:)
      vvar = File->$vvarname$(timeStep,:,:)
    end if
    if (noOfDims .eq. 2) then        ; dims: (x,y)
      uvar = File->$uvarname$(:,:)
      vvar = File->$vvarname$(:,:)
    end if
    if (scaleFactor .ne. 1) then
      uvar = uvar*scaleFactor
      vvar = vvar*scaleFactor
    end if
    if ( useMask ) then
      uvar = mask(uvar,maskVar,0)
      vvar = mask(vvar,maskVar,0)
    end if
    velocity = sqrt(uvar*uvar + vvar*vvar)
    ; set minVar/maxVar for plotting
    if(.not. isvar("minVar")) then
      minVar = min(velocity)
    end if
    if(.not. isvar("maxVar")) then
      maxVar = max(velocity)
    end if
    if ( minVar .gt. maxVar ) then
      print("minVar has to be larger than maxVar")
      exit
    end if
  end if

  ; settings for showgrid
  if (showGrid) then

    print("checkgrid="+checkgrid)

    ; false use of name "cell3" - use attributes instead
    if (checkgrid .eq. "cell3") then
      print("Variable" +varName+ " is located on gridtype " +checkgrid)
      print(" - with showGrid=True plotting of hexagons not yet implemented")
      print(" - routine exits now")
      exit
    end if
 ;  if (plVec)
 ;    print(" - with showGrid=True plotting of Vectors not possible")
 ;    print(" - routine exits now")
 ;    exit
 ;  end if

    boundslonName = File->$lon$@bounds
    boundslatName = File->$lat$@bounds
    boundslon     = File->$boundslonName$ * rad2deg
    boundslat     = File->$boundslatName$ * rad2deg
 ;  use -180 to 180 for mapLLC etc. - error with this correction:
 ;  boundslon = where(boundslon.lt.0, boundslon + 360, boundslon)

    latmax=max(boundslat)
    latmin=min(boundslat)
    lonmax=max(boundslon)
    lonmin=min(boundslon)

    if ( lonmin.lt.0 ) then
    ; longitudes are given in the range [-180,180]
      lonmin = -180.
      lonmax =  180.
    else
    ; longitudes are given in the range [0,360]
      lonmin =    0.
      lonmax =  360.
    end if

    if (DEBUG) then
      print("longitude min/max: " + lonmin + " " + lonmax)  
      print("latitude  min/max: " + latmin + " " + latmax)  
    end if

    ;-----------------------------------------------------------
    ; IMPORTANT: to avoid plotting problem
    ;-----------------------------------------------------------
    ; check longitude

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

        if (DEBUG) then
          print("cell="+icell+" lon="+bl1+", "+bl2+", "+bl3+" changed into " \
                             +" lon="+bn1+", "+bn2+", "+bn3)
        end if

        nfix = nfix +1
      end if
    end do
    print("Longitude(s) of "+nfix+" cells corrected to avoid plotting problem.")

  end if  ;  showgrid

  ;-----------------------------------------------------------
  ; Check output graphics format - rotate eps if necessary
  ;if (oType .eq. "eps") then
  ; the following is not what we want 
  ; ResC@gsnPaperOrientation = "portrait"  ; turns plot on paper
  ;  - we want to rotate the frame orientation including the "paper"
  ;end if

;---------------------------------------------------------------
  if (DEBUG) then
    print("time "+timeVal)
    print("iFile    = "+iFile)
    print("oFile    = "+oFile)
    print("Graphics format is " +oType)
    print("varName  = "+varName)
    print("timeStep = "+timeStep)
    print("mapLLC(lon)   = "+mapLLC(0))
    print("mapLLC(lat)   = "+mapLLC(1))
    print("mapURC(lon)   = "+mapURC(0))
    print("mapURC(lat)   = "+mapURC(1))
    print("minVar   = "+minVar)
    print("maxVar   = "+maxVar)
    if (useMask) then 
      print("maskName   = "+maskName)
    end if
  end if
;---------------------------------------------------------------
;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  ; grid plotting

  ; titles
  if(.not. isvar("lStrg")) then
    dayfrag = timeVal - floor(timeVal)
    day     = doubletointeger(floor(timeVal))
    hour    = doubletointeger(floor(86400*dayfrag/3600))
    minute  = doubletointeger(floor((86400*dayfrag - hour*3600)/60))
    second  = doubletointeger(86400*dayfrag - hour*3600 - minute*60)
    timeStr = str_join((/day,hour,minute,second/),"|")
    if (plotMode.eq."vector") then
      varName = "Velocity"
    end if
    lStrg   = (/varName,"(time "+timeStep+"["+timeStr+"],","pllev "+levIndex+")"/)
  end if
  LeftString   = str_join(lStrg," ")

  if(.not. isvar("rStrg")) then
    rStrg      = (/iFile/)
  end if
  RightString  = str_join(rStrg," ")

  if(.not. isvar("tStrg")) then
    tStrg        = "ICON" 
  end if
  TitleString  = str_join(tStrg," ")
  CenterString = ""

  if(.not. isvar("bStrg")) then
    bStrg      = wcStrt
  end if

  wks = gsn_open_wks(oType,oFile)

  if (showNdc) drawNDCGrid(wks) end if

  gsn_define_colormap(wks,colormap)
  if (showPreview) then
    xwks = gsn_open_wks("x11","x11")
    gsn_define_colormap(xwks,colormap)
  end if

  ResC                      = True

  ResC@gsnMaximize          = False
; ResC@gsnMaximize          = True
; ResC@gsnPaperOrientation  = "landscape"  ; turns plot on paper
; ResC@wkPaperSize          = "A4"
  ResC@gsnFrame             = False
  ResC@gsnSpreadColors      = True

  ResC@cnFillOn             = True   ;  writes labelbar, no contour labels
  ResC@cnLinesOn            = True
;  ResC@cnInfoLabelOn        = False

  FontHeight0               = 0.012
  FontHeight1               = 0.015
  FontHeight2               = 0.017
  FontHeight3               = 0.020

  ResC@tiXAxisFontHeightF   = FontHeight0
  ResC@tiYAxisFontHeightF   = FontHeight0
  ResC@tmXBLabelFontHeightF = FontHeight0
  ResC@tmYLLabelFontHeightF = FontHeight0
  ResC@tmXBLabelJust        = "CenterCenter"
  ResC@gsnStringFontHeightF = FontHeight0



;  ResC@mpFillOn               = True
  if (useMask) then
    ResC@cnMissingValFillColor = -1
  end if

  if (plotMode.eq."scalar") then
    ResC@sfXArray = x      ; triangular mesh
    ResC@sfYArray = y
  end if

  ResC@lbLabelBarOn             = True 
  ResC@pmLabelBarHeightF        = 0.12
  ResC@pmLabelBarWidthF         = 0.80
  ResC@pmLabelBarOrthogonalPosF = 0.15
  ResC@lbLabelAutoStride        = True
  ResC@lbLabelFontHeightF       = FontHeight1     ;   color bar labels

  ResC@gsnLeftString        = LeftString
  ResC@gsnRightString       = RightString
  ResC@gsnCenterString      = CenterString
  ResC@tiMainString         = TitleString


  if (selMode .eq. "manual") then
    ResC@cnLevelSelectionMode = "ManualLevels"
    ResC@cnMinLevelValF       = minVar
    ResC@cnMaxLevelValF       = maxVar
    if (DEBUG) then
      print("Manual Plotlevel setting......")
      print("minVal: " + minVar)
      print("maxVal: " + maxVar)
    end if
    diffLog10                 = log10(abs(maxVar-minVar))
;   if (diffLog10 .lt. 0) then
      ResC@cnLevelSpacingF      = 10^(floor(diffLog10))/contLevs
;   else
;     ResC@cnLevelSpacingF      = 10^(floor(diffLog10))/contLevs
;   end if
  end if
  if (selMode .eq. "halflog") then
    ResC@cnLevelSelectionMode = "ExplicitLevels"
    plotLevels = createLevels(minVar,maxVar,minScale)
    if (DEBUG) print("plotlevels = "+plotLevels) end if
    ResC@cnLevels             = plotLevels
  end if
  ResC@cnFillMode           = "rasterfill"
  ResC@cnRasterSmoothingOn  = True
  ResC@mpGreatCircleLinesOn = True

  if (mapType .eq. "ortho") then
    ResC@mpProjection          = "Orthographic" 
    ResC@mpPerimOn             = False    ; turn off box around plot
    ResC@mpCenterLonF          = 0.    ; specify the plot domain
    ResC@mpCenterLatF          = 0.      ;
    ResC@mpMinLatF             = 0.      ;
    ResC@mpMaxLatF             = 180.      ;

    ResC@mpGridAndLimbOn       = True     ; draw grid lines and limb line
    ResC@mpGridLineDashPattern = 2        ; choose pattern used to draw the grid

  else
    ResC@mpProjection         = "CylindricalEquidistant"
    ResC@mpLimitMode          = "LatLon"
    ResC@mpMinLatF            = mapLLC(1)
    ResC@mpMaxLatF            = mapURC(1)
    ResC@mpMinLonF            = mapLLC(0)
    ResC@mpMaxLonF            = mapURC(0)
    ResC@mpCenterLonF   = ( ResC@mpMaxLonF - ResC@mpMinLonF)/2
  end if

  if (showGrid) then

    plot = gsn_csm_contour_map(wks,var,ResC)
  ; plot = gsn_csm_map(wks,ResC)
    if (showPreview) preview = gsn_csm_contour_map(xwks,var,ResC) end if
    getvalues plot@contour
      "cnLevels"     : levels
      "cnFillColors" : colors
    end getvalues

    flags = new(dimsizes(var),logical,"No_FillValue")
    do i = 0,dimsizes(var) - 1
      flags(i) = where(all(boundslon(i,:) .gt. ResC@mpMaxLonF) .or. \
                       all(boundslon(i,:) .lt. ResC@mpMinLonF) .or. \
                       all(boundslat(i,:) .gt. ResC@mpMaxLatF) .or. \
                       all(boundslat(i,:) .lt. ResC@mpMinLatF), \
                       False, True)
    end do
    if (DEBUG) then
      print(colors)
      print(levels)
      print ("Outside the plot area:           " + dimsizes(ind(flags .eq. False)) + \
             " triangles - not plotted")
    end if

    ; Set up a resource list for the polygons.
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
        if (showPreview) gsn_polygon(xwks,preview,boundslon(vlow(j),:),boundslat(vlow(j),:),pres) end if
      end do
    end if
    if (DEBUG) then
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
        if (showPreview) gsn_polygon(xwks,preview,boundslon(vind(j),:),boundslat(vind(j),:),pres) end if
      end do
      end if
      if (DEBUG) then
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
      if (showPreview) gsn_polygon(xwks,preview,boundslon(vhig(j),:),boundslat(vhig(j),:),pres) end if
    end do
    end if
    if (DEBUG) then
      print ("i=" + i + ", Values   larger than      "+levels(i) +": " \
             + dimsizes(vhig) + " triangles considered, color= " + colors(i+1))
    end if
    delete(vhig)

    ; old try and error - 2011-03
    ; wherehig = where(var .ge. levels(i),True,False)
      ; check index of vhig to decide for plotting the triangle:
      ; print("j="+j+" wherehig="+wherehig(vhig(j))+" vhig="+vhig(j)+" var="+var(vhig))
      ; print("j="+j+" vhig="+vhig(j)+" var="+var(vhig))
;     if (.not. flags(vhig(j))) continue end if
;     if (.not. wherehig(vhig(j))) continue end if

;   wallClockElapseTime(wcStrt, "triangle plot", 0)
    if (showPreview) frame(xwks) end if

  else    ; showGrid
    ; 2011-03-17: error - gsn_vector expects uvar and vvar on lat/lon grid!
    ;  - maybe the attributes could be used to assign the lat/lon of the cell-variables
    ;    on a 2-dim lat/lon variable
    if (plotMode.eq."vector")
     ResC@vcRefMagnitudeF         = 5.0              ; make vectors larger
     ResC@vcRefLengthF            = 8.0            ; ref vec length
     ResC@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
     ResC@vcMinDistanceF          = 0.017            ; thin out vectors

      vc = gsn_csm_vector_scalar_map(wks,uvar,vvar,sqrt(uvar*uvar+vvar*vvar),ResC)
      ;vc = gsn_csm_vector(wks,uvar,vvar,ResC)
      ;vc = gsn_csm_vector_scalar(wks,uvar,vvar,uvar*uvar+vvar*vvar,ResC)
    else  ; shaded contours
      plot  = gsn_csm_contour_map_ce(wks,var,ResC)
      if (showPreview) then
        xplot = gsn_csm_contour_map_ce(xwks,var,ResC)
        frame(xwks)
      end if
    end if
  end if  ; showGrid

  ; base string: time stamp
  ; attention: virtual page changes when lon/lat ratio changes
  txres               = True
  txres@txFontHeightF = 0.011
  gsn_text_ndc(wks,bStrg,.87,.12,txres)

  ; frame overwrites the string?
  frame(wks)
end
