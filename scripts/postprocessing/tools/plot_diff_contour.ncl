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

  varname      = getenv("varname")

  input_file1  = getenv("input_file1")
  input_file2  = getenv("input_file2")

  ; time levels
  istep_start1 = stringtointeger( getenv("time_level_start1") )
  istep_end1   = stringtointeger( getenv("time_level_end1") )
  istep_step1  = stringtointeger( getenv("time_level_step1") )
  istep_start2 = stringtointeger( getenv("time_level_start2") )
  istep_end2   = stringtointeger( getenv("time_level_end2") )
  istep_step2  = stringtointeger( getenv("time_level_step2") )

  ; vertical level
  ilev1        = stringtointeger( getenv("vert_level1") )
  ilev2        = stringtointeger( getenv("vert_level2") )
  if (varname.eq."PS" .or. varname.eq."PHIS") then
    if (ilev1.gt.1) then
      print("WARNING: ignoring ilev1 for a 2-d variable!")
    end if
    if (ilev2.gt.1) then
      print("WARNING: ignoring ilev2 for a 2-d variable!")
    end if
  end if

  ; output options
  output_terminal = getenv("output_terminal")
  output_file = getenv("output_file")

  ; plot options
  proj        = getenv("map_projection")
  colormap    = getenv("colormap")

  ncontours   = stringtointeger( getenv("contour_levels_ncont") )
  if (ncontours.gt.0) then
    contour_start = stringtodouble( getenv("contour_levels_start") )
    contour_end   = stringtodouble( getenv("contour_levels_end") )
  end if

  rel_diff    = stringtointeger( getenv("rel_diff") )
  multiplot   = stringtointeger( getenv("multiplot") )

;---------------------------------------------------------------
; Open data file and get grid information and fields
;---------------------------------------------------------------

  ; Most of the variables must be unique for the difference to make
  ; sense, so they are only read from the first file.
  File  = addfile( input_file1, "r" )
  File2 = addfile( input_file2, "r" )

  ; global attributes of interest
  grid_r      = File@nroot
  grid_b      = File@start_lev
  nlev        = File@nlev
  expname     = File@out_expname
  inpname1    = File@out_expname
  inpname2    = File2@out_expname
  cell        = File@i_cell_type

  ; read the coordinates
  rad2deg     = 45./atan(1.)          ; radians to degrees
  lon         = File->clon *rad2deg   ; longitude [deg]
  lat         = File->clat *rad2deg   ; latitude  [deg]
  times       = File->time            ; time      [days]

  ; Define for each file the three indexes istep_start, istep_end and
  ; nsteps considering that negative start or end indexes indicate the
  ; largest value in the data file.

  ; first file
  if (istep_start1.lt.0) then
    istep_start1 = dimsizes(times)-1
  end if
  if (istep_end1.lt.0) then
    istep_end1 = dimsizes(times)-1
  end if
  if (istep_step1.ne.0) then
    nsteps1 = (istep_end1-istep_start1)/istep_step1 + 1
  else
    nsteps1 = 1
  end if

  ; second file
  if (istep_start2.lt.0) then
    istep_start2 = dimsizes(times)-1
  end if
  if (istep_end2.lt.0) then
    istep_end2 = dimsizes(times)-1
  end if
  if (istep_step2.ne.0) then
    nsteps2 = (istep_end2-istep_start2)/istep_step2 + 1
  else
    nsteps2 = 1
  end if

  ; Total number of steps and consistency checks
  if (nsteps1.eq.nsteps2) then ; OK
    nsteps = nsteps1
  else ; at least one of the two is different from 1
    if(nsteps1.eq.1) then ; OK
      nsteps = nsteps2
    else
      if(nsteps2.eq.1) then ; OK
        nsteps = nsteps1
      else ; WRONG!
        print("ERROR: specified two incompatible ranges!")
        return
      end if
    end if
  end if

  ; read the variables
  isteps1 = new(nsteps,integer)
  isteps2 = new(nsteps,integer)
  var1 = new((/nsteps,dimsizes(lon)/),double)
  var2 = new((/nsteps,dimsizes(lon)/),double)
  do i = 0, nsteps-1
    isteps1(i) = istep_start1 + i*istep_step1
    isteps2(i) = istep_start2 + i*istep_step2
  end do
  do i = 0, nsteps-1
    if (varname.eq."PS" .or. varname.eq."PHIS") then
      var1(i,:) = File ->$varname$(isteps1(i),:)        ; 2-d variable
      var2(i,:) = File2->$varname$(isteps2(i),:)        ; 2-d variable
    else
      var1(i,:) = File ->$varname$(isteps1(i),ilev1-1,:) ; 3-d variable
      var2(i,:) = File2->$varname$(isteps2(i),ilev2-1,:) ; 3-d variable
    end if
  end do

  varlongname = var1@long_name
  if (rel_diff.eq.1) then
    var = (var2-var1)/(max(var2)-min(var2))
    varunit = "rel. err."
  else
    var = var2-var1
    varunit = var1@units
  end if


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
  if(multiplot.ne.1) then ; single plot
    plot = new(1,graphic)
  else
    plot = new(2,graphic)
  end if
  ; make the plot
  gsn_define_colormap(wks,colormap)
  do i = 0, nsteps-1
    if(multiplot.ne.1) then ; single plot
      ResC@gsnCenterString = expname+" step("+isteps2(i)+"-"+isteps1(i)+")"
      plot(0) = gsn_csm_contour_map(wks,var(i,:),ResC)
    else
      ResC@gsnDraw  = False ; do not draw the single plots
      ResC@gsnFrame = False
      
      ResC@lbLabelBarOn = False ; turn off individual colorbars
      ResC@cnLineLabelsOn = False  ; turn off line labels

      ResC@gsnCenterString = inpname1+" step("+isteps1(i)+")"
      plot(0) = gsn_csm_contour_map(wks,var1(i,:),ResC)

      ResC@gsnCenterString = inpname2+" step("+isteps2(i)+")"
      plot(1) = gsn_csm_contour_map(wks,var2(i,:),ResC)

      ResP = True ; resources of the common panel
      ResP@gsnPanelLabelBar = True ; add common colorbar
      gsn_panel(wks,plot,(/2,1/),ResP) ; now draw as one plot
    end if
  end do

end

