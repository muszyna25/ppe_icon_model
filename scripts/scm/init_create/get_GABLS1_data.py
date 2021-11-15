#!/usr/bin/env python
# prepare SCM/LES input and forcing data for the GABLS1
# (marine shallow convection) case in netCDF file
# based on microHH setup of GABLS1 case: https://github.com/microhh/microhh
#
# Ivan Bastak Duran, 2020
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
print ("Beginning of main")

icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_out = icon_dir + data_dir + 'init_SCM_GABLS1.nc'

# specify model levels
hlevels = np.array([0.0,100.,400.,520.,700.,1480.,1500.0,2000.,2100.,3000.,5000.,20000.])
hlevels = hlevels[hlevels.size::-1] # inverse order like model levels ?

# constants
R_d        = 287.06            # dry air gas constant
g          = 9.81              # grav. acceleration
dthetadz   = 0.01
dtdz_str   = 0.04              # stratospheric gradient of potential temperature [K/m]
dtdz_fratm = 0.005             # free atmosphere gradient of potential temperature [K/m]

# variable initialisation
z     = hlevels
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
tke   = np.zeros(np.size(z))
ug    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(z.size):
  # temperature
  if(z[k] <= 100.):
    thl[k] = 265.
  elif(z[k] <= 1500.):
    thl[k] = 265. + dthetadz*(z[k]-100.)
  else:
    thl[k] = 279.+(z[k]-1500.)*dtdz_fratm

  # specific humidity
  qt[:]=0.0
  # u-wind component
  u [:] = 8.

  # TKE
  if(z[k] <= 250.):
    tke[k] = 0.4*(1-z[k]/250.)**3
  else:
    tke[k] =  0.0

  # ug-wind component
  ug[:] = 8.0

  # large scale vertical velocity
  wls[:] = 0.0

  # large scale temperature tendency
  thlls[:] = 0.0

  # large scale moisture tendency
  qtls[:] = 0.0

# normalize profiles to SI
qtls  /= 1000.  # from g/kg to kg/kg
wls   /= 100.   # from cm/s to m/s
thlls /= 86400. # from K/d to K/s
qtls  *= 1.e-8

# surface parameters
p_surf          = 101320 #Pa 
times =np.array([0.0,32400])#s
t_surfs         = np.zeros(times.size)
t_surfs[0]       = 265.#K
t_surfs[1]       = 262.75#K
z0              = 0.1#m
dt_relax        = 7200.0   # s

longitude =  0.0
latitude  =  73.0

#----------------------------------------------------------------

# write data into netcdf file

nlev=hlevels.size
nlevTsoil=9
nlevWsoil=8
nnclass_lu=23
print ("write to netCDF file ", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings for GABLS1 case'
do.createDimension('lev',size=nlev)
do.createDimension('nt', size=times.size)
do.createDimension('levTsoil',size=nlevTsoil)
do.createDimension('levWsoil',size=nlevWsoil)
do.createDimension('nclass_lu',size=nnclass_lu)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units='s since 1969-06-26T00:00:00Z'

#level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=np.arange(nlev)

#soil levels
levTsoil=do.createVariable('levTsoil',np.float64,'levTsoil')
levTsoil.units='level'
levTsoil[:]=np.arange(nlevTsoil)

levWsoil=do.createVariable('levWsoil',np.float64,'levWsoil')
levWsoil.units='level'
levWsoil[:]=np.arange(nlevWsoil)

#atmospheric coordinates
nc_height=do.createVariable('height',np.float64,('nt','lev'))
nc_height.units='meter above surface'

#soil coordinates
nc_depth_T=do.createVariable('depth_T',np.float64,('levTsoil',))
nc_depth_W=do.createVariable('depth_W',np.float64,('levWsoil',))
nc_depth_T.units='depth of T in meter'
nc_depth_W.units='depth of W in meter'

#--- initial profiles

##atmospheric variables
nc_uIN=do.createVariable('uIN',np.float64,('nt','lev'))
nc_uIN.description='zonal wind'
nc_uIN.units='m/s'
nc_vIN=do.createVariable('vIN',np.float64,('nt','lev'))
nc_vIN.description='meridional wind'
nc_vIN.units='m/s'
nc_wIN=do.createVariable('wIN',np.float64,('nt','lev'))
nc_wIN.description='veritcal velocity'
nc_wIN.units='m/s'
nc_thIN=do.createVariable('thIN',np.float64,('nt','lev'))
nc_thIN.description='potential temperature'
nc_thIN.units='K'
nc_qvIN=do.createVariable('qvIN',np.float64,('nt','lev'))
nc_qvIN.description='water vapour specific moisture'
nc_qvIN.units='kg/kg'
nc_qcIN=do.createVariable('qcIN',np.float64,('nt','lev'))
nc_qcIN.description='liquid water specific moisture'
nc_qcIN.units='kg/kg'
nc_pIN=do.createVariable('pIN',np.float64,('nt','lev'))
nc_pIN.description='atmospheric pressure'
nc_pIN.units='Pa'
nc_tkeIN=do.createVariable('tkeIN',np.float64,('nt','lev'))
nc_tkeIN.description='turbulence kinetic energy'
nc_tkeIN.units='m^2/s^2'
###values
nc_height[0,:]=hlevels
nc_time[:]=times
nc_uIN[0,:]=u
nc_uIN[1,:]=u
nc_uIN[:,0]=nc_uIN[:,1]#no gradient of u in stratosphere
nc_vIN[0,:]=0.0
nc_vIN[1,:]=0.0
nc_vIN[:,0]=nc_vIN[:,1]#no gradient of u in stratosphere
nc_wIN[0,:]=0.0
nc_tkeIN[0,:]=tke
nc_tkeIN[1,:]=tke
nc_thIN[0,:]=thl
nc_thIN[1,:]=thl
nc_thIN[:,0]=nc_thIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_qvIN[0,:]=qt
nc_qvIN[1,:]=qt
nc_qvIN[:,0]=nc_qvIN[:,1] #no gradient of moisture in stratosphere
nc_qcIN[:]=0.0
nc_pIN[0,:]=-999.0#should not be used in the simulation
nc_pIN[:,0]=p_surf

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:]=0.0
nc_W_SO[:]=0.0

##forcings
nc_uGEO=do.createVariable('uGEO',np.float64,('nt','lev'))
nc_uGEO.description='zonal geostrophic wind'
nc_uGEO.units='m/s'
nc_vGEO=do.createVariable('vGEO',np.float64,('nt','lev'))
nc_vGEO.description='meridional geostrophic wind'
nc_vGEO.units='m/s'
nc_wLS=do.createVariable('wLS',np.float64,('nt','lev'))
nc_wLS.description='veritcal velocity'
nc_wLS.units='m/s'
nc_dTadv=do.createVariable('dTadv',np.float64,('nt','lev'))
nc_dTadv.description='potential temperature tendency - advection'
nc_dTadv.units='K/s'
nc_dTrad=do.createVariable('dTrad',np.float64,('nt','lev'))
nc_dTrad.description='potential temperature tendency - radiation'
nc_dTrad.units='K/s'
nc_dQVadv=do.createVariable('dQVadv',np.float64,('nt','lev'))
nc_dQVadv.description='potential temperature tendency - advection'
nc_dQVadv.units='kg/kg/s'
nc_dUadv=do.createVariable('dUadv',np.float64,('nt','lev'))
nc_dUadv.description='zonal wind endency - advection'
nc_dUadv.units='m/s/s'
nc_dVadv=do.createVariable('dVadv',np.float64,('nt','lev'))
nc_dVadv.description='meridional wind tendency - advection'
nc_dVadv.units='m/s/s'
nc_height[0,:]=hlevels
nc_height[1,:]=hlevels
###values
nc_uGEO[0,:]=ug
nc_uGEO[1,:]=ug
nc_vGEO[:]=0.0
nc_wLS[:]=0.0
#nc_wLS[:,0]=nc_wLS[:,1]#no gradient in stratosphere
nc_dTadv[:]=0.0
#nc_dTadv[:,0]=nc_dTadv[:,1]#no gradient in stratosphere
nc_dTrad[:]=nc_dTadv[:]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]=0.0
#nc_dQVadv[:,0]=nc_dQVadv[:,1]#no gradient in stratosphere
nc_dUadv[:]=0.0
nc_dVadv[:]=0.0

#--- surface

##pressure
nc_psurf=do.createVariable('psurf',np.float64,('nt'))
nc_psurf.description='Surface atmospheric pressure'
nc_psurf.units='Pa'
##fluxes
nc_sfc_lat_flx=do.createVariable('sfc_lat_flx',np.float64,('nt'))
nc_sfc_lat_flx.description='Latent heat flux'
nc_sfc_lat_flx.units='W/m2'
nc_sfc_sens_flx=do.createVariable('sfc_sens_flx',np.float64,('nt'))
nc_sfc_sens_flx.description='Sensible heat flux'
nc_sfc_sens_flx.units='W/m2'
##temperature and moisture
nc_ts=do.createVariable('ts',np.float64,('nt'))
nc_ts.description='Surface temperature'
nc_ts.units='K'
nc_tg=do.createVariable('tg',np.float64,('nt'))
nc_tg.description='Skin temperature'
nc_tg.units='K'
nc_qvs=do.createVariable('qvs',np.float64,('nt'))
nc_qvs.description='Surface specific vapor humidity'
nc_qvs.units='kg/kg'
## drag coefficients
nc_Cm=do.createVariable('Cm',np.float64,('nt'))
nc_Cm.description='Drag coefficient for heat'
nc_Cm.units=''
nc_Ch=do.createVariable('Ch',np.float64,('nt'))
nc_Ch.description='Drag coefficient for heat'
nc_Ch.units=''
nc_Cq=do.createVariable('Cq',np.float64,('nt'))
nc_Cq.description='Drag coefficient for moisture'
nc_Cq.units=''
## Friction velocity
nc_ustar=do.createVariable('ustar',np.float64,('nt'))
nc_ustar.description='Friction velocity'
nc_ustar.units='m/s'
###values
nc_sfc_lat_flx[:]=0.0
nc_sfc_sens_flx[:]=0.0
nc_psurf[:]=p_surf
nc_ts[:]=t_surfs
nc_tg[:]=t_surfs
nc_qvs[:]=0.0
nc_Ch[:]=0.0
nc_Cq[:]=0.0
nc_Cm[:]=0.0
nc_ustar[:]=0.0

#--- constants

#nudging parameters
nc_dt_relax=do.createVariable('dt_relax',np.float64)
nc_dt_relax.description='relaxation time step for nudging'
nc_dt_relax.units='s'

nc_dt_relax[:]=dt_relax

#external parameters
nc_longitude=do.createVariable('longitude',np.float64)
nc_longitude.description='longitude'
nc_longitude.units='degrees'
nc_latitude=do.createVariable('latitude',np.float64)
nc_latitude.description='latitude'
nc_latitude.units='degrees'
nc_FR_LAND=do.createVariable('FR_LAND',np.float64)
nc_FR_LAND.description='land area fraction'
nc_FR_LAND.units=''
nc_PLCOV_MX=do.createVariable('PLCOV_MX',np.float64)
nc_PLCOV_MX.description='Plant cover maximum due to land use data'
nc_PLCOV_MX.units=''
nc_LAI_MX=do.createVariable('LAI_MX',np.float64)
nc_LAI_MX.description='Leaf Area Index Maximum'
nc_LAI_MX.units=''
nc_ROOTDP=do.createVariable('ROOTDP',np.float64)
nc_ROOTDP.description='Root depth'
nc_ROOTDP.units=''
nc_RSMIN=do.createVariable('RSMIN',np.float64)
nc_RSMIN.description='Minimal stomata resistence'
nc_RSMIN.units=''
nc_SOILTYP=do.createVariable('SOILTYP',np.float64)
nc_SOILTYP.description='soil type derived from FAO Digital Soil Map of the World'
nc_SOILTYP.units=''
nc_Z0=do.createVariable('Z0',np.float64)
nc_Z0.description='Roughness length'
nc_Z0.units=''
nc_EMIS_RAD=do.createVariable('EMIS_RAD',np.float64)
nc_EMIS_RAD.description='longwave surface emissivity'
nc_EMIS_RAD.units=''
nc_TOPO=do.createVariable('TOPO',np.float64)
nc_TOPO.description='geometric height of the earths surface above sea level'
nc_TOPO.units='m'
nc_LU_CLASS_FRACTION=do.createVariable('LU_CLASS_FRACTION',np.float64,('nclass_lu'))
nc_LU_CLASS_FRACTION.description='Fraction of land use classes in target grid element'
nc_LU_CLASS_FRACTION.units=''
nc_LU_CLASS_FRACTION.lctype='GLOBCOVER2009'

nc_longitude[:]=longitude
nc_latitude[:]=latitude
nc_FR_LAND[:]=1.0
nc_PLCOV_MX[:]=0.0
nc_LAI_MX[:]=0.0
nc_ROOTDP[:]=0.0
nc_RSMIN[:]=0.0
nc_SOILTYP[:]=8
nc_Z0[:]=z0
nc_EMIS_RAD[:]=0.98
nc_TOPO[:]=0.0
nc_LU_CLASS_FRACTION[:]=np.zeros([23])
nc_LU_CLASS_FRACTION[5]=1.0#forrest : Closed (>40%) broadleaved deciduous forest (>5m

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()
