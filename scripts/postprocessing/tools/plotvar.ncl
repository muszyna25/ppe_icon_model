
;------------------------------------
; Marco Giorgetta (MPI-M, 2009-05-07)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes color shaded contour plots of a single
; variable for a specified timestep and model level.
;---------------------------------------------------------------

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 begin

;---------------------------------------------------------------
; get environment variables
;---------------------------------------------------------------

  varname     = getenv("icon_varname")

  expfile     = getenv("icon_expfile_1")
  istep       = stringtointeger( getenv("icon_istep_1") )
  ilev        = stringtointeger( getenv("icon_ilev_1") )

  nclwks      = getenv("icon_nclwks")
  plotvarfile = getenv("icon_plotvarfile")

;---------------------------------------------------------------
; open data file and get grid information and fields
;---------------------------------------------------------------

  File        = addfile( expfile, "r" )

; Global attributes of interest
; -----------------------------
  grid_r      = File@nroot
  grid_b      = File@start_lev
  nlev        = File@nlev
  expname     = File@out_expname
  cell        = File@i_cell_type

; Coordinates at cell centers
;----------------------------
  rad2deg     = 45./atan(1.)          ; radians to degrees
  lon         = File->clon *rad2deg   ; longitude [deg]
  lat         = File->clat *rad2deg   ; latitude  [deg]

; Variable at specified timestep and level
;-----------------------------------------
  if (varname.eq."PS" .or. varname.eq."PHIS") then
    var       = File->$varname$(istep,:)        ; 2-d variable
  else
;   level index in file runs from 0 to nlev-1
    ilevm1    = ilev-1
    var       = File->$varname$(istep,ilevm1,:) ; 3-d variable
  end if
;
  varlongname = var@long_name
  varunit     = var@units


;---------------------------------------------------------------
; graphics resources
;---------------------------------------------------------------
  
  ResC                        = True

  ;Contour and shading
  ResC@cnFillOn               = True
  ResC@cnLinesOn              = False

  ;Color bar
  ResC@gsnSpreadColors        = True
  ResC@lbLabelAutoStride      = True

  ;Land sea map
  ResC@mpGeophysicalLineColor = "transparent" ; do not draw land sea map
  ResC@mpFillOn               = False         ; do not fill land sea map

  ;Title string
  ResC@gsnLeftString          = "[1] "+expname+": "+varlongname+" ["+varunit+"], istep="+istep+", ilev="+ilev
  ResC@gsnCenterString        = " "
  ResC@gsnRightString         = "R"+grid_r+"B0"+grid_b+"C"+cell+" L"+nlev

  ;Coordinates of data points
  ResC@sfXArray               = lon
  ResC@sfYArray               = lat

;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks     = gsn_open_wks(nclwks,plotvarfile)
  plot    = new( 1,graphic )
  plot(0) = gsn_csm_contour_map(wks,var,ResC)

end
