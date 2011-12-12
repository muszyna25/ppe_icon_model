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
  nlevs     = cdf_file->var_nlevs(ivar)
  vals      = transpose(cdf_file->values(:, 0:(nlevs-1), ivar, iStation)) ; (time, level, var, station)
  vals@long_name = charactertostring(cdf_file->var_long_name(ivar,:)) \
  +" ("+ charactertostring(cdf_file->var_unit(ivar,:)) + ")"

;----------------------------------------------------------------------
; plot
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
  print("Outputfile '"+str_concat((/oFile,".",oType/))+"' will be created in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens a ps file
  gsn_define_colormap(wks,colormap)  
  res                   = True
  res@tiXAxisString     = "Time Step"              ; x-axis label
  res@tiYAxisString     = "Level"                  ; y-axis label
  res@sfXArray          = cdf_file->time_step(:)
  res@sfYArray          = cdf_file->var_levels(ivar,0:(nlevs-1))
  res@trYReverse        = True                     ; reverses y axis
  res@cnFillOn          = True                     ; turn on color fill
  res@gsnSpreadColors   = True                     ; use full colormap
  res@tiMainFontHeightF = 0.02  

  plot = gsn_contour(wks,vals,res)    ; plot contours
  
end
