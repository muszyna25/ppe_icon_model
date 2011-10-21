;---------------------------------------------------------------
; ICON PLOT LIB
; This library aimes to provide an abstract way to access and plot ICON data
; fields from NetCDF files. It should obay some general requirements:
; * Don-Repeat-Yourself-principle: Every method has one responsebility, which no
;   other method has. Or in other words: Don't just CopyAndPaste.
; * Use as much information from the input files as possible. If something is
;   missing, create a better file or let the user input this peace of
;   information. Examples are coordinates of data variabels and bounds of them.
;   There are attributes for this.
; * naming conventions:
;   # check* methods are procedure which can exit the whole program
;   # get* methods are functions with return, what they descibe with their
;     name. This is similar to NCL itself, e.g. getfilevardims()
;   # set* methods are procedures
;---------------------------------------------------------------
; Authors       : Ralf Mueller (ralf.mueller@zmaw.de)
;                 Stephan Lorenz (stephan.lorenz@zmaw.de)
;-------------------------------------------------------------------------------
; some global defaults
  DEFAULTLON = "clon"
  DEFAULTLAT = "clat"
  RAD2DEG    = 45./atan(1.)
  PLOTMODE   = "scalar"
  ABORTMSG   = "Abort Script."
  NUMLEVS    = 10
;-------------------------------------------------------------------------------
; signum function
undef("sign")
function sign(val)
begin
  if (val.eq.0) then
    return(1)
  else
    return(val/abs(val))
  end if 
end
;-------------------------------------------------------------------------------
; Is the 'time' variable a dimesions in the input dimensions array?
undef("has_time")
function has_time(dims)
begin
  if (any(dims.eq."time")) then
    return(True)
  else
    return(False)
  end if
end
;-------------------------------------------------------------------------------
; Has the input variable the dimension 'time'?
undef("var_has_time")
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
; has a file a certain variable
undef("has_var")
function has_var(filehandle,varname)
begin
  vNames = getfilevarnames(filehandle)
  n      = dimsizes (vNames)
  do i=0,n-1              ; loop thru each variable
    if (vNames(i).eq.varname) return True end if
  end do
  return False
end
;-------------------------------------------------------------------------------
; Read a horizontal field of the variable 'varname' from the given filehandle.
; If it is time dependend or has a vertival axis, use the given values for
; timstep and levelindex for selecting the field.
;
; Currently this limited to 3 dimensions, so that 4D tracer variables are ignored
undef("selIconVar")
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
; see selIconVar(), but for fields on a regular lonlat grid
undef("selRegularField")
function selRegularField(varname,filehandle,timestep,levelindex)
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
; see selRegularField(), but for 3d regular variable
undef("selRegularVar")
function selRegularVar(varname,filehandle,timestep)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
        var = filehandle->$varname$(timestep,:,:,:)
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
        usableVar = False
      end if
    else
      if (noOfDims .eq. 3) then        ; dims: (lev,x,y)
        var = filehandle->$varname$(:,:,:)
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
; scale the given variable with the given factor, if it is NOT 1
undef("scaleVar")
procedure scaleVar(var,scalefactor)
local scalefactor
begin
    if (scalefactor .ne. 1) then
      var = var*scalefactor
    end if
end
;-------------------------------------------------------------------------------
; Return the grid type of a given ICON variable (cell, vertex, or edge).
; This is based on the naming of the coordinates.
undef("getGridType")
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
; Return the coordinates of a given variable in radians
undef("getCoordinates")
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
; Return the coordinates of a given variable in degrees
undef("getLonLats")
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
; Create an ascending array of numbers with a 'halflog' interval, i.e. for
; every decade the numbers 1, 2 and 5 are present.
; Examples:
;   (/1,2,5,10,20,50,100,200,500/) or
;   (/-1,-5e-1,-2e-1,-1e-1,-5e-2,-2e-2,-1e-2,0,1e-2,2e-2,5e-2,1e-1,2e-1,5e-1,1/)
; scaleLimit determines the numer of decades on the positive or negtive axis of the output array
undef("createLevels")
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
; Print all variable names of a given file
undef("printVarNames")
procedure printVarNames(filehandle)
begin
  filevarnames = getfilevarnames(filehandle)
  print("# === variable names in file: "+filevarnames)
end
;---------------------------------------------------------------
; Print information about dimensions and attributes of a given varaiable
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
; Check, if a given variable name is present in a file
undef("checkVarname")
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
; Check, if the dimensions of to variables are equal
undef("checkDimsOfVars")
procedure checkDimsOfVars(varname0,varname1,filehandle)
begin
   print(getfilevardims(filehandle,varname0))
   print(getfilevardims(filehandle,varname1))
  if (str_join(getfilevardims(filehandle,varname0),"") .ne. str_join(getfilevardims(filehandle,varname1),"")) then
    print("Variables '"+varname0+"' and '"+varname1+"' must have the same dimensions!")
    print("ABORT!")
    exit
  end if
end
;---------------------------------------------------------------
; Read a variable from the given or an optional file
undef("getMaskVar")
function getMaskVar(maskvarname,filehandle,hasOtherMaskfile,maskfilename,timestep,levelindex,plotmode)
begin
  if (hasOtherMaskfile) then
    mfilehandle = addfile( maskfilename, "r" )
  else
    mfilehandle = filehandle
  end if
  checkVarname(maskvarname, mfilehandle)
  if (plotmode.eq."scalar" .or. plotmode.eq."overlay") then
    maskvar = selIconVar(maskvarname,mfilehandle,timestep,levelindex)
  else
    if (plotmode.eq."vector") then
      maskvar = selRegularField(maskvarname,mfilehandle,timestep,levelindex)
    else ; section
      maskvar = selRegularVar(maskvarname,mfilehandle,timestep)
    end if
  end if

  return maskvar
end
;---------------------------------------------------------------
; Read the dimensions of a variables from the given ot an optional file
;
; GLOBALS: maskFile
undef("getMaskDim")
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
; Perform a remapping (wich CDO) to a regular grid with a given resolution and return the
; filename of the remapped file. If the file is already present, the name is returned only
undef("remapForVecPlot")
function remapForVecPlot(iFile,resolution,useMask,plotMode,debug)
begin
  if (plotMode.eq."section") then vecVars=(/varName/) end if
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
    variables=maskName+","+str_join(vecVars,",")
    else
    variables=str_join(vecVars,",")
    end if
    if (plotMode.eq."overlay") variables = variables+","+varName end if

    cmd = "cdo -P 8 -remapnn,"+resolution+" -selname,"+variables+" "+iFile+" "+remapFilename
    if debug then
      print(cmd)
    end if
    system(cmd)
  end if
  return remapFilename
end
;---------------------------------------------------------------
; Check, if two files have the same number of timestemps. Exit otherwise
undef("checkRemappedFile")
procedure checkRemappedFile(infilename,remapfilename)
begin
  orgFile   = addfile (infilename , "r")
  remapFile = addfile(remapfilename,"r")
  orgDims   = getvardims(orgFile)
  remapDims = getvardims(remapFile)
  orgTime   = getfilevardimsizes(orgFile,"time")
  remapTime = getfilevardimsizes(remapFile,"time")
  if ( orgTime .eq. remapTime ) then
    print("#=====================================================================================")
    print("Remapped File '"+remapfilename+"' seems to be ok.")
  else
    print("#=====================================================================================")
    print("Remapped File '"+remapfilename+"' seems to be wrong: Has different number of timesteps")
    print("Original file has "+orgTime(0)+" timesteps")
    print("Remapped File has "+remapTime(0)+" timesteps")
    print(ABORTMSG)
    exit
  end if
end
;---------------------------------------------------------------
; Return the bounds of the coordinates of agiven variable
undef("getBoundsOfCoordinates")
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
; (Re)Set end of bounds. Required by the implementation of the grid plot
undef("setBoundsEnds")
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
; Bounds check for the grid plot
undef("checkLongitude")
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
; Set the Vertical Selection mode of the given NCL resource
undef("setLevels")
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
  end if

  if (selmode .eq. "halflog") then
    resource@cnLevelSelectionMode = "ExplicitLevels"
    plotLevels = createLevels(minvar,maxvar,scalelimit)
    if (debug) print("plotlevels = "+plotLevels) end if
    resource@cnLevels             = plotLevels
  end if
  if (debug) then
    print("Manual Plotlevel setting......")
    print("  minVal  : " + minvar)
    print("  maxVal  : " + maxvar)
    print("  numLevs: " + numlevs)
    print("  scalelimit: " + scalelimit)
  end if
end
;---------------------------------------------------------------
; Create a default NCL resource
undef("setDefaultResource")
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
; Compute the left shift for string with font hight 0.01, which should be
; displayed in the lower left corner of the plot with gsn_test_ndc()
undef("getLshiftForNDCString")
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
; Display a small string in the lower left corner of the plot
undef("setBaseString")
procedure setBaseString(workstation,plot,resource,basestring)
begin
  resource@txFontHeightF = 0.01 ; if this is changed, please adjust the implementation of getLshiftForNDCString()
  gsn_text_ndc(workstation,basestring,getLshiftForNDCString(basestring),.01,resource)
end
;---------------------------------------------------------------
; see setBaseString, but use an optional global argument
;
; GLOBALS: bStrg
undef("setAutomaticBaseString")
procedure setAutomaticBaseString(workstation,plotmode)
begin
  if(.not. isvar("bStrg")) then
    bStrg      = ""
  end if
  basestr     = str_join(bStrg," ")
  ; base string: plot file name + basestring + time stamp
  ; attention: virtual page changes when lon/lat ratio changes
  BaseString = "Prgr icon_ocean.ncl: " + " " + basestr + " " + systemfunc("date")
  res               = True
  res@txFontHeightF = 0.011
  if (plotmode.eq."section") res@txPosYF = 0.01 end if
  gsn_text_ndc(workstation,BaseString,getLshiftForNDCString(BaseString),.12,res)
end
;---------------------------------------------------------------
; Set the Captions of a plot (left,right,center,title)
undef("setPlotCaptions")
procedure setPlotCaptions(resource,leftstring,rightstring,centerstring,titlestring)
begin
  resource@gsnLeftString   = leftstring
  resource@gsnRightString  = rightstring
  resource@gsnCenterString = centerstring
  resource@tiMainString    = titlestring
end
;---------------------------------------------------------------
; Set the capitons of a plot automatically
;
; GLOBALS: lStrg, rStrg, tStrg
undef("setAutomaticPlotCaptions")
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
; Determine the type of the map
undef("setMapType")
procedure setMapType(resource,maptype,centerlon,centerlat,satdist)
begin
  if (maptype .eq. "ortho") then
    resource@mpProjection              = "Orthographic"
    resource@mpPerimOn                 = False;     turn off box around plot
    resource@mpGridAndLimbOn           = True;      draw grid lines and limb line
    resource@mpGridLineDashPattern     = 2;         choose pattern used to draw the grid
    resource@mpCenterLonF              = centerlon
    resource@mpCenterLatF              = centerlat
    resource@vcRefLengthF              = 0.01
    resource@vcRefAnnoArrowLineColor   = "black"; change ref vector color
    resource@vcRefAnnoArrowUseVecColor = False;  do not use vec color for ref
    resource@vcRefAnnoOn               = False;        turns off ref vector annotation
    resource@vcGlyphStyle              = "CurlyVector"; turn on curly vectors
    resource@vcMonoLineArrowColor      = True;        multiple colors desired
    resource@vcLineArrowColor          = "blue"; change vector color
    resource@vcVectorDrawOrder         = "PostDraw"; draw vectors last


    resource@mpGridAndLimbOn        = True               ; turn on lat/lon grid lines
    resource@mpGridMaskMode         = "MaskNotOcean"     ; don't draw over land or
    ; inland water bodies

;    resource@mpLandFillColor        = "tan"
;    resource@mpOceanFillColor       = "LightBlue"
;    resource@mpInlandWaterFillColor = "LightBlue"
  end if

  if (maptype .eq. "lonlat")
    resource@mpProjection         = "CylindricalEquidistant"
    resource@mpLimitMode          = "LatLon"
  end if

  if (maptype .eq. "NHps")
    resource@mpProjection        = "Stereographic"
    resource@mpRelativeCenterLon = True
    resource@mpCenterLonF        = 0
    resource@mpRelativeCenterLat = True
    resource@mpCenterLatF        = 90.
    resource@mpLimitMode         = "Corners"
  end if

  if (maptype .eq. "SHps")
    resource@mpProjection        = "Stereographic"
    resource@mpRelativeCenterLon = True
    resource@mpCenterLonF        = 0
    resource@mpRelativeCenterLat = True
    resource@mpCenterLatF        = -90.
    resource@mpLimitMode         = "Corners"
  end if
  if (maptype .eq. "sat")
    resource@mpProjection     = "Satellite"
    resource@mpMinLatF        = resource@mpMinLatF + 10.
    resource@mpMaxLatF        = resource@mpMaxLatF - 10.
    resource@mpMinLonF        = resource@mpMinLonF + 10.
    resource@mpMaxLonF        = resource@mpMaxLonF - 10.
    resource@mpCenterLonF     = centerlon
    resource@mpCenterLatF     = centerlat
    resource@mpSatelliteDistF = satdist
  end if
end
;---------------------------------------------------------------
; Determine the corners of the map to plot
; mapllc = map-lower-left-corner
; mapurl = map-upper-right-corner
; format: (/lon,lat/)
undef("selMapCut")
procedure selMapCut(resource,mapllc,mapurc)
begin
  ; bounds of plotting area (not used for Orthographic)
  resource@mpMinLatF    = mapllc(1)
  resource@mpMaxLatF    = mapurc(1)
  resource@mpMinLonF    = mapllc(0)
  resource@mpMaxLonF    = mapurc(0)

  ; center of view (for Orthographic and Sattelite only)
;  resource@mpCenterLonF = (resource@mpMaxLonF+resource@mpMinLonF)/2.0
;  resource@mpCenterLatF = (resource@mpMaxLatF+resource@mpMinLatF)/2.0

  ; bounds of plotting area for stereographic plots
  resource@mpLeftCornerLonF = mapllc(0)
  resource@mpLeftCornerLatF = mapllc(1)
  resource@mpRightCornerLonF = mapurc(0)
  resource@mpRightCornerLatF = mapurc(1)
end
;---------------------------------------------------------------
; Set the visibility of the map
undef("setMapVisibility")
procedure setMapVisibility(resource,mapline)
begin
  if (mapline) then
    resource@mpGeophysicalLineColor = "foreground"
  else
    resource@mpGeophysicalLineColor = "transparent"
  end if
end
;---------------------------------------------------------------
; exit if minvar is larger than maxvar
undef("checkMinMaxVar")
procedure checkMinMaxVar(minvar,maxvar)
begin
  if ( minvar .gt. maxvar ) then
    print("minVar has to be larger than maxVar")
    exit
  end if
end
;---------------------------------------------------------------
; set the color of masked locations
undef("setMaskColor")
procedure setMaskColor(wks,resource)
begin
  newcolor = NhlNewColor(wks,0.85,0.85,0.85)   ; add gray to colormap
  ;resource@gsnSpreadColorStart  = 1
  resource@gsnSpreadColorEnd     = 19  ; last color number not known?
  resource@cnMissingValFillColor = 20  ; last color is used for missing values
end
;---------------------------------------------------------------
; set the coordinates for an ICON data field
undef("setCoordinates")
procedure setCoordinates(resource,lons,lats)
begin
  resource@sfXArray = lons
  resource@sfYArray = lats
end
;---------------------------------------------------------------
; Print out some information about the used map
undef("showMapInfo")
procedure showMapInfo(resource,maptype,mapline)
begin
  print("mapType: "+maptype)
  print("mapLine: "+mapline)
  print("CentLon/Lat="+resource@mpCenterLonF+"; "+resource@mpCenterLatF)
  print("Min/MaxLatF="+resource@mpMinLatF+"; "+resource@mpMaxLatF)
  print("Min/MaxLonF="+resource@mpMinLonF+"; "+resource@mpMaxLonF)
end
;---------------------------------------------------------------
; setting for vertical cross sections based on remapped icon data
undef("setSection")
function setSection(secLC,secRC,npoints,var)
begin
  leftlat  = secLC(1)
  rightlat = secRC(1)

  leftlon  = secLC(0)
  rightlon = secRC(0)

  dist     = gc_latlon(leftlat,leftlon,rightlat,rightlon,npoints,2)
  trans    = linint2_points(var&lon,var&lat,var,True,dist@gclon,dist@gclat,2)

  return trans
end
;---------------------------------------------------------------
; Create flags array the the grid plot
undef("getFlags")
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
; Set the default vector plot resource
undef("setDefaultVectorPlot")
procedure setDefaultVectorPlot(resource, refmagnitude, reflength, style, mindistance)
begin
  resource@vcRefMagnitudeF         = refmagnitude ; make vectors larger
  resource@vcRefLengthF            = reflength    ; reference vector length
  resource@vcGlyphStyle            = style        ; turn on curly vectors
  resource@vcMinDistanceF          = mindistance  ; thin out vectors
end
;---------------------------------------------------------------
; Return the Colors of a given plot
undef("getColorsFromPlot")
function getColorsFromPlot(plot)
begin
  getvalues plot@contour
    "cnFillColors" : colors
  end getvalues
  return colors
end
;---------------------------------------------------------------
; Return the levels of a given plot
undef("getLevelsFromPlot")
function getLevelsFromPlot(plot)
begin
  getvalues plot@contour
    "cnLevels"     : levels
  end getvalues
  return levels
end
;---------------------------------------------------------------
; Print out some information about the grid plot
undef("printGridPlotInfo")
procedure printGridPlotInfo(colors,levels,flags)
begin
  print(colors)
  print(levels)
  print ("Outside the plot area:           " + dimsizes(ind(flags .eq. False)) + \
      " triangles - not plotted")
end
;---------------------------------------------------------------
; Create a resource for the grid plot
undef("setupGridResource")
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
; Plot the ICON grid of a variable to the plot.
undef("plotGrid")
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
;---------------------------------------------------------------
; read in mask variable
undef("readMaskVar")
function readMaskVar()
begin
end
