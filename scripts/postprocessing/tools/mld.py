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
iconDir  = '/pool/data/ICON/tools'
if os.environ.has_key('ICONDIR'):
  iconDir = os.environ.get('ICONDIR')
iconPlot = iconDir+'/icon_plot.ncl'
iconLib  = iconDir
hostname = os.popen('hostname').read().strip()

#=============================================================================
# CDO setup
cdo.forceOutput = not 'FORCE' in os.environ
cdo.debug       = not 'DEBUG' in os.environ
cdoOptions      = '-r'

#=============================================================================
# output meta data
file      = io.open("partab","w")
tagString = u"&PARAMETER\n CODE=18\n NAME=mixed_layer_depth\n STANDARD_NAME=mixed_layer_depth\n LONG_NAME='Mixed layer depth'\n UNITS='m'\n/"
file.write(tagString)
file.close()

# ONLY USE MARCH FOR NORTHERN HEMISPHERE, SEPTEMBER for SOUTHERN, SEPTEMBER for SOUTHERN
#unless cdo.showmon(:input => ifile)[0].split.include?('3')
#  warn "Could not find march in input data!"
# exit(1)
#end unless ENV['CHECK'].nil?

# check is T and S are present
names = cdo.showname(input=ifile)[0].split()


def checkVars(allVarNames):
  # basic setup for t and s
  conf               = {}
  conf['name_t']     = 't'
  conf['name_s']     = 's'
  conf['use_t']      = True
  conf['use_s']      = True
  conf['use_tAcc']   = False
  conf['use_sAcc']   = False
  conf['maskName']   = 'wet_c'
  conf['useMask']    = False
  conf['ignoreTime'] = os.environ.has_key('IGNORE_TIME')

  keys = [x for x in conf.keys() if x[0:4] == 'name']
  varnames2Check = []
  for k in keys:
    varnames2Check.append(k[5:])
# print(varnames2Check)
  for name in varnames2Check:
    if name not in allVarNames:
      conf['use_'+name] = False
      # check for accumulatied values
      if name+'_acc' not in allVarNames:
        conf['use_'+name+'Acc'] = False
        print('could not find temperature (t,t_acc)')
      else:
        conf['use_'+name+'Acc'] = True
        conf['name_'+name] = name+'_acc'
    else:
      print('variable '+name+' not in file')

  if conf['maskName'] in allVarNames:
    conf['useMask'] = True

  return conf

check  = checkVars(names)
tsSelection, maskName, selTime = '','',ifile


inputArg = ''
if check['use_tAcc'] and check['use_sAcc']:
  tsSelection = '-chname,t_acc,t,s_acc,s -selname,t_acc,s_acc'
if check['use_t'] and check['use_s']:
  tsSelection = '-selname,t,s'
if not check['ignoreTime']:
  selTime = '-selmon,3,9 '
if check['useMask']:
  maskVar  = "-selname,%s -seltimestep,1"%(check['maskName'])
  inputArg = "-adisit -setcode,-1 -div %s %s %s %s %s"%(tsSelection,selTime,ifile,maskVar,ifile)
else:
  inputArg = "-adisit -setcode,-1 %s %s %s "%(tsSelection,selTime,ifile)

# select T and S, set code to -1 to be ignored by rhopot/adisit
cdo.rhopot(0,
           input   = inputArg,
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
otype    = '-oType=png'
otype    = ''
title    = "-tStrg=%s"%(ifile)
for t in range(0,ntime):
  t = str(t)
  os.system("nclsh %s -altLibDir=%s -iFile=%s -varName=mixed_layer_depth -oFile=mld_%s -timeStep=%s %s %s %s %s\n"%(iconPlot,iconLib,mld,t,t,select,colormap,otype,title))
