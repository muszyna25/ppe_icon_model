
;------------------------------------
; Marco Giorgetta (MPI-M, 2009-05-07)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes color shaded contour plots of a single
; variable for specified timesteps and model levels.
; The script plots 3 panels: 
; [1] variable from source 1 (file 1, level 1, time step 1)
; [2] variable from source 2 (file 2, level 2, time step 2)
; [3] difference [2]-[1] 
;---------------------------------------------------------------

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 begin

;---------------------------------------------------------------
; get environment variables
;---------------------------------------------------------------

  varname     = getenv("icon_varname")

; source 1
; --------
  expfile_1   = getenv("icon_expfile_1")
  istep_1     = stringtointeger( getenv("icon_istep_1") )
  ilev_1      = stringtointeger( getenv("icon_ilev_1") )

; source 2
; --------
  expfile_2   = getenv("icon_expfile_2")
  istep_2     = stringtointeger( getenv("icon_istep_2") )
  ilev_2      = stringtointeger( getenv("icon_ilev_2") )

  nclwks      = getenv("icon_nclwks")
  plotdifffile= getenv("icon_plotdifffile")

;---------------------------------------------------------------
; open data file and get grid information and fields
;---------------------------------------------------------------

  File_1  = addfile( expfile_1, "r" )
  File_2  = addfile( expfile_2, "r" )

; Global attributes of interest
; -----------------------------
  grid_r      = File_1@nroot
  grid_b      = File_1@start_lev
  nlev        = File_1@nlev
  expname_1     = File_1@out_expname
  expname_2     = File_2@out_expname
  cell        = File_1@i_cell_type

; Coordinates at cell centers
;----------------------------
  rad2deg = 45./atan(1.)          ; radians to degrees
  lon     = File_1->clon *rad2deg ; longitude [deg]
  lat     = File_1->clat *rad2deg ; latitude  [deg]

; Variable at specified timesteps and levels
;-------------------------------------------
  if (varname.eq."PS" .or. varname.eq."PHIS") then  ; 2-d variable
    var_1     = File_1->$varname$(istep_1,:)
    var_2     = File_2->$varname$(istep_2,:)
  else                                                ; 3-d variable
;   level index in file runs from 0 to nlev-1
    ilevm1_1  = ilev_1-1
    ilevm1_2  = ilev_2-1
    var_1     = File_1->$varname$(istep_1,ilevm1_1,:)
    var_2     = File_2->$varname$(istep_2,ilevm1_2,:)
  end if

  varlongname = var_1@long_name
  varunit     = var_1@units

  var_3 = var_2-var_1

;---------------------------------------------------------------
; graphics resources
;---------------------------------------------------------------
  
; Panel 1
; -------
  ResC_p1                        = True

  ;Contour and shading
  ResC_p1@cnFillOn               = True
  ResC_p1@cnLinesOn              = False

  ;Color bar
  ResC_p1@gsnSpreadColors        = True
  ResC_p1@lbLabelAutoStride      = True

  ;Land sea map
  ResC_p1@mpGeophysicalLineColor = "transparent" ; do not draw land sea map
  ResC_p1@mpFillOn               = False         ; do not fill land sea map

  ;Title string
  ResC_p1@gsnLeftString          = "[1] "+expname_1+": "+varlongname+" ["+varunit+"], istep="+istep_1+", ilev="+ilev_1
  ResC_p1@gsnCenterString        = " "
  ResC_p1@gsnRightString         = "R"+grid_r+"B0"+grid_b+"C"+cell+" L"+nlev

  ;Coordinates of data points
  ResC_p1@sfXArray               = lon
  ResC_p1@sfYArray               = lat

  ; Do not draw plots immediately, but only after "gsn_panel"
  ResC_p1@gsnDraw                = False
  ResC_p1@gsnFrame               = False


; Panel 2
; -------
  ResC_p2                        = True
  ResC_p2                        = ResC_p1 ; copy attributes from ResC_p1

  ;Title string
  ResC_p2@gsnLeftString          = "[2] "+expname_2+": "+varlongname+" ["+varunit+"], istep="+istep_2+", ilev="+ilev_2


; Panel 3
; -------
  ResC_p3                        = True
  ResC_p3                        = ResC_p1 ; copy attributes from ResC_p1

  ;Title string
  ResC_p3@gsnLeftString          = "[3]=[2]-[1]: " + varlongname + " [" + varunit + "]"


;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks     = gsn_open_wks(nclwks,plotdifffile)
  plot    = new( 3,graphic )
  plot(0) = gsn_csm_contour_map(wks,var_1,ResC_p1)
  plot(1) = gsn_csm_contour_map(wks,var_2,ResC_p2)
  plot(2) = gsn_csm_contour_map(wks,var_3,ResC_p3)


;------------------------------------------
; panel
;------------------------------------------

  ResP = True
  gsn_panel (wks, plot, (/3,1/), ResP)

end
