;----------------------------------------------------------------------
; mtgrm_plot.ncl 
;----------------------------------------------------------------------
; Meteogram plot script for vertical variables.
;
; Basic Usage:
;   ncl -n mtgrm_plot_sfc.ncl 'iFile="path/file.nc"' iStation=1 'varName="SHFL"' 'expnum="${expnum}"'
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
  nvars = dimsizes(cdf_file->sfcvar_name)
  do ivar=0,(nvars(0)-1)
    iname = charactertostring(cdf_file->sfcvar_name(ivar,:))
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

  xstride      = 6  ; subdivisions of X-Axis labels (small ticks)

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(iFile,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable
  ivar      = get_var_index(cdf_file, varName)
  vals      = cdf_file->sfcvalues(0:, ivar, iStation)   ; (time, var, station)
  time      = cdf_file->time_step(0:)

;----------------------------------------------------------------------
; plot setup
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
  print("Outputfile '"+str_concat((/oFile,".000001.",oType/))+"' will be created in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens a ps file
  gsn_define_colormap(wks,colormap)  
  date = convert_date(charactertostring(cdf_file->date(0,:)))
  yyyymmdd = sprinti("%0.2i", date@year)+sprinti("%0.2i", date@month)+sprinti("%0.2i", date@day)

  res                   = True
  res@gsnCenterString   = charactertostring(cdf_file->sfcvar_long_name(ivar,:))
  res@gsnCenterStringOrthogonalPosF = 0.06
  res@gsnLeftString     = expnum
  res@gsnRightString    = charactertostring(cdf_file->station_name(iStation,:))
  res@tiXAxisString     = "Time  (from "+yyyymmdd+")"                                  ; x-axis label
  res@tiYAxisString     = "["+ charactertostring(cdf_file->sfcvar_unit(ivar,:)) + "]"  ; y-axis label
  res@gsnCenterStringFontHeightF = 0.02
  res@gsnLeftStringFontHeightF   = 0.018
  res@gsnRightStringFontHeightF  = 0.018
  res@tiYAxisFontHeightF         = 0.018
  res@tiXAxisFontHeightF         = 0.018
  res@tmXBLabelFontHeightF       = 0.014
  res@vpXF                 = 0.15                  ; viewport left margin (def: 0.2)
  res@vpWidthF             = 0.75                  ; The Width of the plot box
; res@vpHeightF            = 0.58                  ; The height of the plot box
  res@tmXTOn               = False                 ; turn off the top tickmarks
  res@xyLineThicknesses    = 5                     ; increase line thickness
  res@xyLineColor          = "blue"                ; set line color

;----------------------------------------------------------------------
; create values and date labels for time axis explicitly
;----------------------------------------------------------------------

  ntimesteps = dimsizes(cdf_file->time_step(:))
  expl_labels = new(ntimesteps(0), string)
  expl_values = new(ntimesteps(0), float)

  count = 0
  do i=0,(ntimesteps(0)-1)  
    if (mod(i,xstride) .eq. 0) then
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
; line plot
;----------------------------------------------------------------------

  plot = gsn_csm_xy(wks,time,vals,res)             ; plot line
  
end
