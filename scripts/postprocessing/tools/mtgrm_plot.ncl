;----------------------------------------------------------------------
; mtgrm_plot.ncl
;----------------------------------------------------------------------
; Example plot of data in meteogram NetCDF.
; Note: We assume here that we deal with a 3D variable, main levels.
; Invocation:
;   ncl -n mtgrm_plot.ncl
;----------------------------------------------------------------------
; 12/2011 : F. Prill, DWD
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Utility function-----------------------------------------------------
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
; Constants
;----------------------------------------------------------------------
  DataFileName = "METEOGRAM_patch001.nc"
  variable     = "P" 
  istation     = 0
  itime        = 0
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(DataFileName,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable
  ivar      = get_var_index(cdf_file, variable)
  nlevs     = cdf_file->var_nlevs(ivar)
  vals      = transpose(cdf_file->values(:, :, ivar, istation)) ; (time, level, var, station)
  vals@long_name = charactertostring(cdf_file->var_long_name(ivar,:)) \
  +" ("+ charactertostring(cdf_file->var_unit(ivar,:)) + ")"
  
; plot header --------------------------------------------------------
  print( "Meteogram '" + charactertostring(cdf_file->station_name(istation,:))+"'" )
  
  wks  = gsn_open_wks("ps","mtgrm_example") ; Opens a ps file
  gsn_define_colormap(wks,"BlWhRe")  
  res                   = True
  res@tiXAxisString     = "time step"              ; x-axis label
  res@tiYAxisString     = "level"                  ; y-axis label
  res@sfXArray          = cdf_file->time_step(:)
  res@sfYArray          = cdf_file->var_levels(ivar,0:(nlevs-1))
  res@trYReverse        = True                     ; reverses y axis
  res@cnFillOn          = True                     ; turn on color fill

  res@gsnSpreadColors   = True                     ; use full colormap"
  res@tiMainFontHeightF = 0.02  

  plot = gsn_contour(wks,vals,res)    ; plot contours
  
end
