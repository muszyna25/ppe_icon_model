;----------------------------------------------------------------------
; mtgrm_plot.ncl 
;----------------------------------------------------------------------
; Meteogram plot script for vertical variables.
;
; Basic Usage:
;   ncl -n mtgrm_plot.ncl 'iFile="path/file.nc"' iStation=1 'varName="T"'
;
; Required Parameters:
; iFile         : input file (e.g. "METEOGRAM_patch001.nc")
; iStation      : station number (1,2,3,...)
; varName       : name of the variable to (contour) plot
;
; Optional Parameter:
; oFile         : plot file without extension (set by the output type: oType)
; oType         : output graphics format (ps, eps, png, default: eps)
; colormap      : string for predefined colormaps of the plot (e.g. 'colormap="BlAqGrYeOrReVi200"')
;
; Info:
; * There is a meteogram plot script available publically, which might be used to add features.
;   http://www.ncl.ucar.edu/Applications/meteo.shtml
; * Defaults: Note that the default escape character is colon ":", which is
;   changed to "~" in our .hluresfile.  Also the default font is set to helvetica.
;
;----------------------------------------------------------------------
; 12/2011 : F. Prill,   DWD (florian.prill@dwd.de)
; 12/2011 : M. Koehler, DWD (martin.koehler@dwd.de)
;----------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;----------------------------------------------------------------------
; Function to find variable in output list
;----------------------------------------------------------------------

function get_var_index (cdf_file, zname : string)
local iname, ivar, nvars
begin
  nvars = dimsizes(cdf_file->var_name)
  do ivar=0,(nvars(0)-1)
    iname = charactertostring(cdf_file->var_name(ivar,:))
    if(zname .eq. iname) then
      return(ivar)
    end if
  end do
  return -1
end


;----------------------------------------------------------------------
; Function to get substring from string
;----------------------------------------------------------------------

function substring (instr[1]:string, first[1]:numeric, last[1]:numeric)
local instr, first, last, main, p2
begin
   main = stringtochar (instr)   
   if (last .ge. first) then
      p2 = last
   else
      p2 = dimsizes (main) - 2
   end if
   return (chartostring (main(first:p2)))
end


;----------------------------------------------------------------------
; Function to retrieve date from iso8601 string
;    Input:  string in the format of "yymmddThhmmssZ"
;    Output: year, month, day, hour, minute, second
;----------------------------------------------------------------------

function convert_date (zstr : string)
begin
  date = True
  date@year   = stringtointeger(substring(zstr, 0, 3))
  date@month  = stringtointeger(substring(zstr, 4, 5))
  date@day    = stringtointeger(substring(zstr, 6, 7))
  date@hour   = stringtointeger(substring(zstr, 9,10))
  date@minute = stringtointeger(substring(zstr,11,12))
  date@second = stringtointeger(substring(zstr,13,14))
  return date
end


begin

;----------------------------------------------------------------------
; argument list
;----------------------------------------------------------------------

  if(.not. isvar("oType")) oType = "eps" end if
  if(.not. isvar("oFile")) then
    ext      = get_file_suffix(iFile,0)
    oFile    = ext@fBase
  end if
  if (.not. isvar("colormap"))    colormap    = "BlueDarkRed18" end if

  iStation     = iStation - 1
  itime        = 0

  stride       = 6  ; subdivisions of X-Axis labels
  ystride      = 10 ; subdivisions of right Y-Axis labels

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(iFile,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable
  ivar      = get_var_index(cdf_file, varName)
  nlevs     = cdf_file->var_nlevs(ivar)
  vals      = transpose(cdf_file->values(:, 0:(nlevs-1), ivar, iStation)) ; (time, level, var, station)
  vals@long_name = charactertostring(cdf_file->var_long_name(ivar,:))

;----------------------------------------------------------------------
; plot setup
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
  print("Outputfile '"+str_concat((/oFile,".",oType/))+"' will be created in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens a ps file
  gsn_define_colormap(wks,colormap)  
  date = convert_date(charactertostring(cdf_file->date(0,:)))
  yyyymmdd = sprinti("%0.2i", date@year)+sprinti("%0.2i", date@month)+sprinti("%0.2i", date@day)

  res                   = True
  res@gsnRightString    = charactertostring(cdf_file->station_name(iStation,:))
  res@tiXAxisString     = "Time  (from "+yyyymmdd+")" ; x-axis label
  res@tiYAxisString     = "Model Level"               ; y-axis label
  res@sfXArray          = cdf_file->time_step(:)
  res@sfYArray          = cdf_file->var_levels(ivar,0:(nlevs-1))
  res@trYReverse        = True                        ; reverses y axis
  res@cnFillOn          = True                        ; turn on color fill
  res@gsnSpreadColors   = True                        ; use full colormap
  res@gsnLeftStringFontHeightF  = 0.02
  res@gsnRightStringFontHeightF = 0.017
  res@tiYAxisFontHeightF   = 0.018
  res@tiXAxisFontHeightF   = 0.018

;----------------------------------------------------------------------
; contour and label bar resources
;----------------------------------------------------------------------

  res@lbTitleOn            = True                  ; turn on title
  res@lbTitleString        = "["+ charactertostring(cdf_file->var_unit(ivar,:)) + "]"
  res@lbTitleFontHeightF   = .015                  ; label title font size
 ;res@lbTitlePosition      = "Right"               ; label position
  res@lbOrientation        =   "vertical"          ; vertical label bar
  res@pmLabelBarWidthF     = 0.24                  ; horizontal shift
  res@lbBoxMinorExtentF    = 0.15                  ; horizontal width
  res@lbLabelFontHeightF   = 0.015                 ; label font size
 ;res@cnLinesOn            = False                 ; optional lines
 ;res@cnFillMode           = "RasterFill"          ; optional raster
 ;res@cnLevelSelectionMode = "ExplicitLevels"
 ;res@gsnMaximize          = True                  ; use full page
  res@vpXF                 = 0.11                  ; viewport left margin (def: 0.2)
  res@vpWidthF             = 0.69                  ; viewport width (def: 0.6)

;----------------------------------------------------------------------
; create values and date labels for time axis explicitly
;----------------------------------------------------------------------

  ntimesteps = dimsizes(cdf_file->time_step(:))
  expl_labels = new(ntimesteps(0), string)
  expl_values = new(ntimesteps(0), float)

  count = 0
  do i=0,(ntimesteps(0)-1)  
    if (mod(i,stride) .eq. 0) then
      date = convert_date(charactertostring(cdf_file->date(i,:)))
      date_str = ""
      if ((i .eq. 0) .or. (.not. (date@day .eq. last_date@day))) then
        date_str = sprinti("%0.2i", date@month) \
        +          "-" + sprinti("%0.2d", date@day) + " "
      end if
      expl_values(count) = cdf_file->time_step(i)
      expl_labels(count) = date_str + sprinti("%0.2d", date@hour) \
      +                   ":" + sprinti("%0.2d", date@minute)
      count = count + 1
      last_date = date
    end if
  end do

  res@tmXBMode             = "Explicit"            ; Define own tick mark labels.
  res@tmXBMinorOn          = True
  res@tmXBMinorValues      = cdf_file->time_step(:)
  res@tmXBMinorThicknessF  = 0.1
  res@tmXBValues           = expl_values(0:(count-1))  ; location of explicit labels
  res@tmXBLabels           = expl_labels(0:(count-1))
  res@tmXBLabelAngleF      = 90.             ; change label angle
  res@tmXBLabelFontHeightF = .01             ; Font size
  res@tmXBLabelJust        = "CenterRight"
  ; remove mirrored tick marks (top, right):
  res@tmXTOn               = False
  res@tmYROn               = False

;----------------------------------------------------------------------
; create an independent axis with level heights
;----------------------------------------------------------------------

  expl_labels2 = new(nlevs, string)
  do i=0,(nlevs-1)
    expl_labels2 = sprintf("%3.0f m", cdf_file->heights(0:(nlevs-1), ivar, iStation))
  end do
  xy3 = create "xy3" xyPlotClass wks
    "trYReverse"              : "True"
    "trYMaxF"                 :  cdf_file->var_levels(ivar, nlevs-1)
    "trYMinF"                 :  cdf_file->var_levels(ivar, 0)
    "tmYROn"                  : "True"
    "tmYLOn"                  : "False"
    "tmYUseLeft"              : "False"
    "tmYRLabelsOn"            : "True"
    "tmYRMajorLengthF"        : .01
    "tmYRMajorOutwardLengthF" : .01
    "tmXBOn"                  : "False"
    "tmXTOn"                  : "False"
    "tmYLOn"                  : "False"
    "tmYROn"                  : "True"
    "tmYRMode"                : "Explicit"
    "tmYRValues"              : cdf_file->var_levels(ivar,0:(nlevs-1):ystride)
    "tmYRLabels"              : expl_labels2(0:(nlevs-1):ystride)
    "tmYRLabelJust"           : "CenterRight"
    "tiYAxisSide"             : "Right"
    "tmYRLabelFontHeightF"    : 0.01
    "tiXAxisFontHeightF"      : 0.01
    "tiYAxisFontHeightF"      : 0.03
    "tmYRMinorOn"             : "False"
    "tmYLMinorOn"             : "False"
    "tmXBMinorOn"             : "False"
    "tmXTBorderOn"            : "False"
    "tmXBBorderOn"            : "False"
    "tmYLBorderOn"            : "False"
  end create
  draw(xy3)

;----------------------------------------------------------------------
; contour plot
;----------------------------------------------------------------------

  plot = gsn_csm_contour(wks,vals,res)
  
end
