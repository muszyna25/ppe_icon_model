;---------------------------------------------------------------
; This script makes contour/vector plots of general ICON data files
; For scalar variables the underlying grid can be used instead of automatic
; contour lines. Both modes are capable of masking the input files before
; plotting, see the command line options below for details.
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
; oFile         : plot file wihtout extension (it set by the output type: oType)
; varName       : name of the variable to (contour) plot
; vecVars       : space separated string of the 2 vector components u and v to draw a vector plot
;
; Optional Parameter:
; oType         : output graphics format (ps, eps, png, default: eps)
; resolution    : resolution string whish is used for remapping the icon data
;                 to regular grid (default: r90x45)
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
; scaleLimit    : Limits the number of decades for levels in selMode-halflog
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
; showGrid      : display polygon lines with filled colors instead of contour plot
; showNcd       : display NDC Grid to find Normalized Device Coordinates on the plot
;---------------------------------------------------------------
; altLibDir     : Alternative directory for loading icon_plot_lib.ncl
;---------------------------------------------------------------
; Authors       : Ralf Mueller (ralf.mueller@zmaw.de)
;                 Stephan Lorenz (stephan.lorenz@zmaw.de)
;---------------------------------------------------------------

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

 
 if ( .not. isvar("altLibDir") )  then
   print("load default library from /pool/data/ICON/tools")
   loadscript("/pool/data/ICON/tools/icon_plot_lib.ncl")
 else
   print("Load library: "+altLibDir+"/icon_plot_lib.ncl")
   loadscript(altLibDir +"/icon_plot_lib.ncl")
 end if

;===============================================================
;----- MAIN PROGRAM --------------------------------------------
;- options handling + checking
 begin
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

  gsn_define_colormap(wks,colormap)

  if (showPreview) then
    xwks = gsn_open_wks("x11","x11")
    gsn_define_colormap(xwks,colormap)
  end if

  if (showNdc)
    drawNDCGrid(wks)
    if (showPreview) drawNDCGrid(xwks) end if
  end if

  ResC = True

  setDefaultResource(ResC)

  if (maxView) ResC@gsnMaximize = True end if

  if (useMask) setMaskColor(wks,ResC) end if

  if (plotMode.eq."scalar") setCoordinates(ResC,x,y) end if

  setAutomaticPlotCaptions(ResC,plotMode,varName,File,iFile,timeStep,levIndex)

  setAutomaticBaseString(wks)
  if (showPreview) setAutomaticBaseString(xwks) end if

  setLevels(selMode,ResC,minVar,maxVar,scaleLimit,numLevs,DEBUG)

  setMapType(ResC,mapType)

  selMapCut(ResC,mapLLC,mapURC)

  setSateliteProj(ResC,projSat)

  setMapVisibility(ResC,mapLine)

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
    if (DEBUG) print("Gridtype is "+getGridType(var)) end if
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
