----------------------------------------------------------------------------
Notes on MARS4ICON / IFS2ICON usage
----------------------------------------------------------------------------
Initial version: F.Prill, DWD (2011-12-16)

Document version:
  $LastChangedDate$
  $Rev: 7638 $ -- $Author$

----------------------------------------------------------------------------

Table of contents:

1. GRIB_API
2. MARS Retrieval
3. IFS2ICON script
3.1 Precomputation of weights
3.2 IFS2ICON usage


1. GRIB API
-----------

Necessary:  grib_api Version 1.9.9
Check, eg., with "grib_info"; cf. settings in ".profile".

Important: Use grib_api *release* version
           This is necessary for setting "CI" missing values to 0 !
             (. ${GRIB_API}/bin/grib_api-setenv -r -V1.9.9)

Table of GRIB keys:
  http://www.ecmwf.int/publications/manuals/d/gribapi/param
see also
  firefox --no-remote /usr/local/pkg/grib_api/1.8.1/html/GRIB_PARAMETER_html/index.htm )

The IFS2ICON scripts handles GRIB1 input fields by their table codes
and GRIB2 input fields by their short names.

Get GRIB info of retrieved file:
  grib_ls ifs_oper_T255_2011010100.grb

Get GRIB version (GRIB1/GRIB2) of a given file:
  grib_get  -w count=1 -p editionNumber ifs_oper_T1279_2011080100.grb



2. MARS Retrieval
-----------------

Script for retrieval of IFS data (via MARS, ECMWF): See ICON svn repository,

~/trunk/icon-dev/scripts/preprocessing/mars4icon_smi_32r3+ 

Run as:
  mars4icon_smi_32r3+ 2011010100 0 255 1/to/91 2>&1 | tee log.mars4icon.txt

Note: Requires valid ECMWF certificate (command "eccert" on lxe)

Note: At the end of this script, volumetric soil water layer fields (SWVL1-4) are
      renamed into "soil moisture index" (SMIL1); coded as GRIB indices 80-83. This
      is essential because values are converted without altering the grib index,
      according to the procedure described in

      http://www.ecmwf.int/products/changes/soil_hydrology_cy32r3/#compute_smi



3. IFS2ICON script
------------------

3.1 Precomputation of weights
-----------------------------

CDO command for precomputation of weights: "cdo gencon" (conservative remapping),
as descibed in

Jones, P. W. First- and Second-Order Conservative Remapping Schemes for Grids in
 Spherical Coordinates Mon. Wea. Rev., 1999, 127, 2204-2210

MUST USE /e/uhome/extrmuel/local/bin/cdo  !!!
(Check version, at least 1.5.3, with ".../cdo -V"

Execute 
 >    /e/uhome/extrmuel/local/bin/cdo -P 8 gencon,/e/uhome/gzaengl/icon-1.0.6_RESTRUCT/experiments/GRF_R2B6N8_EUR2/iconR2B06_DOM01.nc ifs_T255_2011010100_new.grb cell_weight_ifs255_R2B06.nc
before first application of IFS2ICON script!

Perform this on lxt, see also Redmine Features #1712 (performance issues).



3.2 IFS2ICON usage
------------------

Show info for Ruby script version:
     fprill@lxt:~/IFS2ICON> ifs2icon.rb -v

Prepare config file,
e.g. "config/ifs2icon_hori_smi_255_to_R2B6_DOM01_speedtest.conf", especially
- INPUTFILE
- OUTPUTFILE
- GRIDFILE
- CELLWEIGHTFILE
- OUTPUTVARS
  Note: 2D variable "CI" (sea ice) is necessary for ICON but uses different set of
  weights (this take some more hours...)

Note: OPTIONS FOR HYDROSTATIC OUTPUT - unused but must be set anyway
Note: VERBOSE and DEBUG will be set to FALSE for later (operational) use

Then apply IFS2ICON:
     ./ifs2icon.rb -c <config/config_file.conf>

Note: The ifs2icon script distinguishes between GRIB edition 1 and 2 by
      parsing the file extension, "*.grb" or "*.grb2" !

Note: Usually the variable "CI" (sea ice) takes MUCH longer for interpolation.
      Reason: Weights are computed "on the fly" because this 2D field contains
      missing values corresponding to a land-sea mask..
      Due to a fix by M. Koehler (2011-12-16), this bottleneck has been 
      removed, because GRIB fields from MARS are downloaded without mask, and 
      missing values are replaced by "0".

