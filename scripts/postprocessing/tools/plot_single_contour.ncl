;------------------------------------
; NCL script for plotting a single scalar variable. 
; 
; WARNING: this script should be called through the driver
; "plot_scalar_and_vector", it should not be launched directly from
; the user. See plot_scalar_and_vector -h for details.
;
; Marco Restelli (MPI-M, 2009-03-31)
;------------------------------------
; Marco Giorgetta (MPI-M, 2009-03-31)
; - read in global attribute and print also experiment name
;   and grid resolution in title bar

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 begin

;---------------------------------------------------------------
; Get environment variables
;---------------------------------------------------------------

  varname     = getenv("varname")

  input_file  = getenv("input_file")

  ; time levels
  istep_start = stringtointeger( getenv("time_level_start") )
  istep_end   = stringtointeger( getenv("time_level_end") )
  istep_step  = stringtointeger( getenv("time_level_step") )

  ; vertical level
  ilev        = stringtointeger( getenv("vert_level") )
  if (varname.eq."PS" .or. varname.eq."PHIS") then
    if (ilev.gt.1) then
      print("WARNING: ignoring ilev for a 2-d variable!")
    end if
  end if

  proj        = getenv("map_projection")
  colormap    = getenv("colormap")

  ncontours   = stringtointeger( getenv("contour_levels_ncont") )
  if (ncontours.gt.0) then
    contour_start = stringtodouble( getenv("contour_levels_start") )
    contour_end   = stringtodouble( getenv("contour_levels_end") )
  end if

  output_terminal = getenv("output_terminal")
  output_file = getenv("output_file")

;---------------------------------------------------------------
; Open data file and get grid information and fields
;---------------------------------------------------------------

  File = addfile( input_file, "r" )

  ; global attributes of interest
  grid_r      = File@nroot
  grid_b      = File@start_lev
  nlev        = File@nlev
  expname     = File@out_expname
  cell        = File@i_cell_type

  ; read the coordinates
  rad2deg     = 45./atan(1.)          ; radians to degrees
  lon         = File->clon *rad2deg   ; longitude [deg]
  lat         = File->clat *rad2deg   ; latitude  [deg]
  times       = File->time            ; time      [days]
  ; negative start or end indexes indicate the largest value
  if (istep_start.lt.0) then
    istep_start = dimsizes(times)-1
  end if
  if (istep_end.lt.0) then
    istep_end = dimsizes(times)-1
  end if
  nsteps = (istep_end-istep_start)/istep_step + 1

  ; read the variable
  isteps = new(nsteps,integer)
  var = new((/nsteps,dimsizes(lon)/),double)
  do i = 0, nsteps-1
    isteps(i) = istep_start + i*istep_step
  end do
  do i = 0, nsteps-1
    if (varname.eq."PS" .or. varname.eq."PHIS") then
      var(i,:) = File->$varname$(isteps(i),:)        ; 2-d variable
    else
      var(i,:) = File->$varname$(isteps(i),ilev-1,:) ; 3-d variable
    end if
  end do

  varlongname = var@long_name
  varunit     = var@units

;---------------------------------------------------------------
; Setup resources
;---------------------------------------------------------------

  ResC = True

  ; coordinates of data points
  ResC@sfXArray = lon
  ResC@sfYArray = lat

  ResC@cnFillOn  = True  ; use colors
  ResC@cnLinesOn = False ; do not mark contour lines

  ResC@lbLabelAutoStride = True ; fix colorbar

  ResC@gsnSpreadColors = True ; use the whole color bar
  if (ncontours.gt.0) then
    contour_levs = new(ncontours,double)
    if (ncontours.gt.1) then
      contour_interval = (contour_end-contour_start)/(ncontours-1)
      do i = 0, ncontours-2
        contour_levs(i) = contour_start + i*contour_interval
      end do
    end if
    ; make sure the last element is exact
    contour_levs(ncontours-1) = contour_end
    ResC@cnLevelSelectionMode = "ExplicitLevels" 
    ResC@cnLevels = contour_levs
  end if

  ResC@mpProjection = proj  ; map projection
  ResC@mpPerimOn    = False ; do not plot the bounding box

  ; do not plot the land sea map
  ResC@mpGeophysicalLineColor = "transparent" ; do not draw land sea map
  ResC@mpFillOn               = False         ; do not fill land sea map

  ResC@gsnLeftString          = varlongname+" ["+varunit+"]"
  ResC@gsnRightString         = "R"+grid_r+"B0"+grid_b+"C"+cell+" L"+nlev

;---------------------------------------------------------------
; Make the plot
;---------------------------------------------------------------

  ; create the workspace and 
  wks  = gsn_open_wks(output_terminal,output_file)
  ; make the plot
  gsn_define_colormap(wks,colormap)
  do i = 0, nsteps-1
    ResC@gsnCenterString = expname+" step"+isteps(i)
    plot = gsn_csm_contour_map(wks,var(i,:),ResC)
  end do

end
