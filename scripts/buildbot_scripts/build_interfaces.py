#! /usr/bin/env python
# -*- coding: utf-8 -*-
#==============================================================================
import os
from model_paths import *

paths=model_paths()
#-----------------------------------------------------------------------


def make_all_binaries(configure_flags):
  os.chdir(paths.basePath)
  status = os.system("./configure "+configure_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = os.system("./build_command")
  if not status == 0:
    print("Build failed")
    return status
  os.chdir(paths.thisPath)
  return 0
    
def make_ocean_binaries(configure_flags):
  ocean_only_flags=" --disable-atmo --disable-jsbach  --with-yac=no"
  os.chdir(paths.basePath)
  ocean_folder = paths.basePath+"/ocean"
  status = os.system("scripts/building/get_ocean "+paths.basePath+" "+ocean_folder)
  if not status == 0:
    print("get_ocean failed")
    return status
  os.chdir(ocean_folder)  
  status = os.system("./configure "+configure_flags+ocean_only_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = os.system("./build_command")
  if not status == 0:
    print("Build failed")
    return status

  # get the ocean binaries and setup-info  
  os.chdir(paths.basePath)
  status = os.system("cp -r "+ocean_folder+"/build .")
  status = os.system("cp "+ocean_folder+"/config/set-up.info config")
  
  return 0


# if succesful returns a list of the runscripts
#  otherwise returns the status
def make_runscript(experimentPathName, runflags):
  os.chdir(paths.basePath)
  # seperate the the input path from the experiment name
  experimentPath, experimentName = paths.getPathAndName(experimentPathName)
  # separate prefix and main name
  experimentPrefixName = experimentName.split(".",2)
  #print(experimentPath+" "+experimentPrefixName[0]+" "+experimentPrefixName[1])
  outscript=experimentName+".run"
  make_runscript_command=paths.basePath+"/config/make_target_runscript "
  inoutFiles="in_script="+experimentPathName+" in_script=exec.iconrun out_script="+outscript+" EXPNAME="+experimentPrefixName[1]+" "
  status = os.system(make_runscript_command+inoutFiles+runflags)
  if not status == 0:
    print("make_runscripts failed")
  os.chdir(paths.thisPath)
  return status, outscript

#-----------------------------------------------------------------------

