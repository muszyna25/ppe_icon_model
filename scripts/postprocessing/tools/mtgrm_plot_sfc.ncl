;----------------------------------------------------------------------
; mtgrm_plot.ncl 
;----------------------------------------------------------------------
; Meteogram plot script for vertical variables.
;
; Basic Usage:
;   ncl -n mtgrm_plot.ncl 'iFile="path/file.nc"' iStation=1 'varName="SHFL"'
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

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(iFile,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable
  ivar      = get_var_index(cdf_file, varName)
  vals      = cdf_file->sfcvalues(0:, ivar, iStation) ; (time, level, var, station)
  time      = cdf_file->time_step(0:)
  ntime     = dimsizes(time)
 ;print(ntime)
 ;print(time)
 ;print(vals)

;----------------------------------------------------------------------
; plot
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
  print("Outputfile '"+str_concat((/oFile,".",oType/))+"' will be created in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens a ps file
  gsn_define_colormap(wks,colormap)  
  res                   = True
  res@tiMainString      = charactertostring(cdf_file->sfcvar_long_name(ivar,:))
  res@tiYAxisString     = "["+ charactertostring(cdf_file->sfcvar_unit(ivar,:)) + "]"  ; y-axis label
  res@tiXAxisString     = "Time"                                                       ; x-axis label

  res@tmXBMode          = "Explicit"               ; Define own tick mark labels.
  res@tmXBValues        = time(0:ntime-1:24)       ; location of explicit labels
  res@tmXBMinorValues   = time(0:ntime-1)          ; location of explicit labels
  res@tmXBLabels        = charactertostring(cdf_file->date(1:ntime-1:24,:))
  res@tmLabelAutoStride = True

  res@tiMainFontHeightF    = 0.02 
  res@tiYAxisFontHeightF   = 0.018
  res@tiXAxisFontHeightF   = 0.018
  res@tmXBLabelFontHeightF = 0.014
  res@vpWidthF          = 0.70   ; The Width of the plot box
  res@vpHeightF         = 0.50   ; The height of the plot box
; res@trXMaxF           = 72.0   ; max value on x-axis
; res@trYMaxF           = max(tempht)   ; max value on y-axis
; res@trYMinF           = min(tempht)   ; min value on y-axis
  
; res@tmXBMode          = "Explicit" ; Define own tick mark labels.
; res@tmXBValues        = taus       ; location of explicit labels
; res@tmXBLabels        = (/"12z", "", "18z",  "", "Apr29", \
;                    "", "06z", "",  "12z", "",   \
;                    "18z", "", "Apr30","", "06z",   \
;                    "", "12z", "",  "18z", "",   \
;                    "May01", "", "06z",  "", "12z"/)
  res@tmXTOn            = False         ; turn off the top tickmarks
  res@xyLineThicknesses = 2          ; increase line thickness
  res@xyLineColor       =  "blue"        ; set line color

  plot = gsn_csm_xy(wks,time,vals,res)   ; plot line
  
end
