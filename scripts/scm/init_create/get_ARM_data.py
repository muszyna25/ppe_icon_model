#!/usr/bin/env python
# prepare SCM/LES input and forcing data for the BOMEX
# (marine shallow convection) case in netCDF file
# based on microHH setup of ARM case: https://github.com/microhh/microhh
#
# Ivan Bastak Duran, 2019
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
print ("Beginning of main")

icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_out = icon_dir + data_dir + 'init_SCM_ARM.nc'

# specify model levels
hlevels = np.array([0.0,50.,350.,650.,700.,1000.,1300.,2500.,5500.,10000,20000.])
hlevels = hlevels[hlevels.size::-1] # inverse order like model levels ?

# constants
R_d        = 287.06            # dry air gas constant
g          = 9.81              # grav. acceleration
dtdz_str   = 0.04              # stratospheric gradient of potential temperature [K/m]
dtdz_fratm = 0.005             # free atmosphere gradient of potential temperature [K/m]

# variable initialisation
z     = hlevels
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
v     = np.zeros(np.size(z))
tke   = np.zeros(np.size(z))
ug    = np.zeros(np.size(z))
vg    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(z.size):
  # temperature
  if(z[k] <= 50.):
    thl[k] = 299.0  + (z[k]     )*(301.5 -299.0 )/(50.)
    qt[k] = 15.20  + (z[k]     )*(15.17 -15.20 )/(50.)
  elif(z[k] <=  350.):
    thl[k] = 301.5  + (z[k]-  50.)*(302.5 -301.5 )/(350.-50.)
    qt[k] = 15.17  + (z[k]-  50.)*(14.98 -15.17 )/(350.-50.)
  elif(z[k] <=  650.):
    thl[k] = 302.5  + (z[k]- 350.)*(303.53-302.5 )/(650.-350.)
    qt[k] = 14.98  + (z[k]- 350.)*(14.80 -14.98 )/(650.-350.)
  elif(z[k] <=  700.):
    thl[k] = 303.53 + (z[k]- 650.)*(303.7 -303.53)/(700.-650.)
    qt[k] = 14.80  + (z[k]- 650.)*(14.70 -14.80 )/(700.-650.)
  elif(z[k] <= 1300.):
    thl[k] = 303.7  + (z[k]- 700.)*(307.13-303.7 )/(1300.-700.)
    qt[k] = 14.70  + (z[k]- 700.)*( 13.50-14.80 )/(1300.-700.)
  elif(z[k] <= 2500.):
    thl[k] = 307.13 + (z[k]-1300.)*(314.0 -307.13)/(2500.-1300.)
    qt[k] = 13.50  + (z[k]-1300.)*( 3.00 - 13.50)/(2500.-1300.)
  elif(z[k] <= 5500.):
    thl[k] = 314.0  + (z[k]-2500.)*(343.2 -314.0 )/(5500.-2500.)
    qt[k] =  3.00
  else:
    thl[k] = 343.2 + (z[k]-5500.)*dtdz_fratm

  # u-wind component
  u[:] = 10.

  # ug-wind component
  ug[k] = 10.

# normalize profiles to SI
qt /= 1000.  # g to kg

# set the time series
t  = np.array([  0.,   4.,  6.5,  7.5,  10., 12.5, 14.5])
tout  = np.arange(0,15.0,0.5)

H  = np.interp(tout,t,np.array([-30.,  90., 140., 140., 100., -10.,  -10]))
LE = np.interp(tout,t,np.array([  5., 250., 450., 500., 420., 180.,    0]))

advthl = np.array([ 0.   , 0.  ,  0.  , -0.08, -0.16, -0.16])
radthl = np.array([-0.125, 0.  ,  0.  ,  0.  ,  0.   , -0.1])
advqt  = np.array([ 0.08 , 0.02, -0.04, -0.10, -0.16, -0.30])

tls    = np.array([  0.,   3.,  6.,  9.,  12., 14.5])
thlls_in  = np.zeros((t.size, z.size))
qtls_in   = np.zeros((t.size, z.size))
thlls  = np.zeros((tout.size, z.size))
qtls   = np.zeros((tout.size, z.size))

# large scale forcings
for n in range(tls.size):
  tendthl = advthl[n] + radthl[n]
  tendqt  = advqt[n]
  for k in range(z.size):
    # temperature
    if(z[k] <= 1000.):
      thlls_in[n,k] = tendthl
      qtls_in[n,k] = tendqt
    else:
      thlls_in[n,k] = tendthl - (z[k]-1000.)*(tendthl)/(5500.-1000.)
      qtls_in[n,k] = tendqt - (z[k]-1000.)*(tendqt)/(5500.-1000.)
    thlls[:,k]  = np.interp(tout,t,thlls_in[:,k])
    qtls[:,k]  = np.interp(tout,t,qtls_in[:,k])


#unit conversion
t           *= 3600.   # h to s
tout        *= 3600.   # h to s
tls         *= 3600.   # h to s
thlls       /= 3600.   # h to s
thlls_in    /= 3600.   # h to s
qtls        /= 3600.   # h to s
qtls        /= 1000.   # g to kg
qtls_in     /= 3600.   # h to s
qtls_in     /= 1000.   # g to kg


# surface parameters
ustar        = 0.0     # m/s    --no slip for ARM
sfc_lat_flx  = LE      # W/m2
sfc_sens_flx = H       # W/m2
p_surf       = 97000.  # surface pressure [Pa]
z0           = 0.035   # m
dt_relax     = 7200.0  # s

longitude = 0.0
latitude  = 36.0

#----------------------------------------------------------------

# write data into netcdf file

nlev=hlevels.size
nlevTsoil=9
nlevWsoil=8
nnclass_lu=23
print ("write to netCDF file ", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings for ARM case'
do.createDimension('lev',size=nlev)
do.createDimension('nt', tout.size)
do.createDimension('levTsoil',size=nlevTsoil)
do.createDimension('levWsoil',size=nlevWsoil)
do.createDimension('nclass_lu',size=nnclass_lu)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units='s since 1997-06-21T00:00:00Z'

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
nc_time[:]=tout
for it in  range(0,tout.size): #in this case all are the same
  nc_height[it,:]=hlevels
  nc_uIN[it,:]=u
  nc_uIN[it,0]=nc_uIN[it,1]#no gradient of u in stratosphere
  nc_vIN[it,:]=v
  nc_vIN[it,0]=nc_vIN[it,1]#no gradient of u in stratosphere
  nc_wIN[it,:]=0.0
  nc_tkeIN[it,:]=tke
  nc_thIN[it,:]=thl
  nc_thIN[it,0]=nc_thIN[it,1]+dtdz_str*(nc_height[it,0]-nc_height[it,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
  nc_qvIN[it,:]=qt
  nc_qvIN[it,0]=nc_qvIN[it,1] #no gradient of moisture in stratosphere
  nc_qcIN[it,:]=0.0
  nc_pIN[it,:]=-999.0#should not be used in the simulation
  nc_pIN[it,0]=p_surf

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
###values
for it in  range(0,tout.size):
 nc_uGEO[it,:]=ug
 nc_vGEO[it,:]=vg
 nc_wLS[it,:]=wls
 nc_wLS[it,0]=nc_wLS[it,1]#no gradient in stratosphere

nc_dTadv[:,:]=thlls
nc_dTadv[:,0]=nc_dTadv[it,1]#no gradient in stratosphere
nc_dTrad[:,:]=nc_dTadv[:,:]*0.0 # setting radiation tendency to 0
nc_dQVadv[:,:]=qtls
nc_dQVadv[:,0]=nc_dQVadv[it,1]#no gradient in stratosphere
nc_dUadv[:,:]=0.0
nc_dVadv[:,:]=0.0

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
nc_sfc_lat_flx[:]=sfc_lat_flx
nc_sfc_sens_flx[:]=sfc_sens_flx
nc_psurf[:]=p_surf
nc_ts[:]=thl[-1]
nc_tg[:]=thl[-1]
nc_qvs[:]=0.0
nc_Ch[:]=0.0
nc_Cq[:]=0.0
nc_Cm[:]=0.0
nc_ustar[:]=ustar

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
