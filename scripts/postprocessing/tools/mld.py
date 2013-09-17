#!/usr/bin/env python
from cdo import *
from jobqueue import *
import os,sys,io

def dbg(msg):
  if 'DEBUG' in os.environ:
    print msg
#=============================================================================
# input setup, names for intermediate files
# icon input
if len(sys.argv) == 1:
  print "No input file given"
  sys.exit(1)
else:
  ifile=sys.argv[1]

rhopot       = 'rhopot.nc'
rhopot_delta = 'rhopot_deltaToSurface.nc'
mld          = 'icon_mld.nc'
q            = JobQueue(8)
cdo          = Cdo()
#cdo.setCdo('cdo-dev')

#=============================================================================
# plotting setup
iconplot = '/pool/data/ICON/tools/icon_plot.ncl'
iconlib  = '/pool/data/ICON/tools'
hostname = os.popen('hostname').read().strip()
plot     = iconplot + ' ' + iconlib

#=============================================================================
# CDO setup
cdo.forceOutput = not 'FORCE' in os.environ
cdo.debug       = not 'DEBUG' in os.environ
cdoOptions      = ''

#=============================================================================
# output meta data
file      = io.open("partab","w")
tagString = u"&PARAMETER\n CODE=18\n NAME=mixed_layer_depth\n STANDARD_NAME=mixed_layer_depth\n LONG_NAME='Mixed layer depth'\n UNITS='m'\n/"
file.write(tagString)
file.close()

# select T and S, set code to -1 to be ignored by rhopot/adisit
# ONLY USE MARCH FOR NORTHERN HEMISPHERE
#unless cdo.showmon(:input => ifile)[0].split.include?('3')
#  warn "Could not find march in input data!"
# exit(1)
#end unless ENV['CHECK'].nil?
cdo.rhopot(0,
           input   = "-adisit -setcode,-1 -div -selname,T,S -selmon,3 %s -selname,wet_c %s"%(ifile,ifile),
           output  = rhopot,
           options = cdoOptions)

# substracto the surface value
cdo.sub(input  = "%s -sellevidx,1 %s"%(rhopot,rhopot),
        output = rhopot_delta)
# compute the depth if the iso surface for a value of 0.125
cdo.setpartab('partab',
              input  = "-isosurface,0.125 %s"%(rhopot_delta),
              output = mld)
#=============================================================================
# plot each timestep
ntime    = int(cdo.ntime(input = mld)[0])
# select north atlantic
select   = '-mapLLC=-60,30 -mapURC=30,85'
colormap = '-colormap=testcmap'
for t in range(0,ntime):
  t = str(t)
  #print("nclsh %s -iFile=%s -varName=mixed_layer_depth -oFile=mld_%s -isIcon -timeStep=%s %s %s \n"%(plot,mld,t,t,select,colormap))
  os.system("nclsh %s -iFile=%s -varName=mixed_layer_depth -oFile=mld_%s -isIcon -timeStep=%s %s %s \n"%(plot,mld,t,t,select,colormap))
