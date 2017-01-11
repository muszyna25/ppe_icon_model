#! /usr/bin/env python
# -*- coding: utf-8 -*-
#==============================================================================
import os,subprocess
from model_paths import paths

#-----------------------------------------------------------------------
# some help method to call commands and return stdout and stderr
def runCommand(cmd):
# proc = subprocess.Popen(' '.join(cmd),
#                         shell  = True,
#                         stderr = subprocess.PIPE,
#                         stdout = subprocess.PIPE)
#
# retvals = proc.communicate()
# stdout  = retvals[0].decode("utf-8")
# stderr  = retvals[1].decode("utf-8")

  print("Command:"+cmd)
# print("StdOut:"+stdout)
# print("StdErr:"+stderr)

  status = subprocess.call(cmd,shell=True)
  return status


#-----------------------------------------------------------------------

def make_binaries_interface(configure_flags, builder_flags):
  if "Ocean" in builder_flags:
    return make_ocean_binaries(configure_flags)
  #elif "AES" in builder_flags:
    #return make_aes_binaries(configure_flags)
  elif "NWP" in builder_flags:
    return make_nwp_binaries(configure_flags)
  elif not "Inactive" in builder_flags:
    return make_all_binaries(configure_flags)
  return 0
  

def make_all_binaries(configure_flags):
  os.chdir(paths.basePath)
  status = runCommand("./configure "+configure_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = runCommand("./build_command")
  if not status == 0:
    print("Build failed")
    return status
  set_account()
  os.chdir(paths.thisPath)
  return 0
    
def make_ocean_binaries(configure_flags):
  ocean_flags=" --disable-atmo --disable-jsbach --with-yac=no --with-flags=ocean"
  os.chdir(paths.basePath)
  ocean_folder = paths.basePath+"/ocean_build"
  status = runCommand("scripts/building/get_ocean "+paths.basePath+" "+ocean_folder)
  if not status == 0:
    print("get_ocean failed")
    return status
  os.chdir(ocean_folder)  
  status = runCommand("./configure "+configure_flags+ocean_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = runCommand("./build_command")
  if not status == 0:
    print("Build failed")
    return status
  # get the ocean binaries and setup-info  
  os.chdir(paths.basePath)
  status = runCommand("cp -r "+ocean_folder+"/build .")
  status = runCommand("cp "+ocean_folder+"/config/set-up.info config")
  set_account()
  
  os.chdir(paths.thisPath)  
  return 0

def make_aes_binaries(configure_flags):
  aes_flags=" --disable-ocean --with-yac=no"
  os.chdir(paths.basePath)
  aes_folder = paths.basePath+"/aes_build"
  status = runCommand("scripts/building/get_atmo "+paths.basePath+" "+aes_folder)
  if not status == 0:
    print("get_aes failed")
    return status
  os.chdir(aes_folder)
  status = runCommand("./configure "+configure_flags+aes_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = runCommand("./build_command")
  if not status == 0:
    print("Build failed")
    return status
  # get the aes binaries and setup-info
  os.chdir(paths.basePath)
  status = runCommand("cp -r "+aes_folder+"/build .")
  status = runCommand("cp "+aes_folder+"/config/set-up.info config")
  set_account()

  os.chdir(paths.thisPath)
  return 0


def make_nwp_binaries(configure_flags):
  aes_flags=" --disable-ocean --disable-jsbach --disable-testbed --with-yac=no"
  os.chdir(paths.basePath)
  aes_folder = paths.basePath+"/nwp_build"
  status = runCommand("scripts/building/get_atmo "+paths.basePath+" "+aes_folder)
  if not status == 0:
    print("get_aes failed")
    return status
  os.chdir(aes_folder)
  status = runCommand("./configure "+configure_flags+aes_flags)
  if not status == 0:
    print("Configure failed")
    return status
  status = runCommand("./build_command")
  if not status == 0:
    print("Build failed")
    return status
  # get the aes binaries and setup-info
  os.chdir(paths.basePath)
  status = runCommand("cp -r "+aes_folder+"/build .")
  status = runCommand("cp "+aes_folder+"/config/set-up.info config")
  set_account()

  os.chdir(paths.thisPath)
  return 0

def set_account():
  setup_file = open("./config/set-up.info", 'a')
  setup_file.write("use_account_no=mh0287\n")
  setup_file.write("use_notification=never\n")
  setup_file.close()
  
# if succesful returns a list of the runscripts
#  otherwise returns the status
def make_runscript(experimentPathName, runflags):
  if not paths.thisExperimentExists(experimentPathName):
    print("Warning:experiment "+experimentPathName+" does not exist.")
    return 1
  #print("make runscript:"+experimentPathName+" with flags:"+runflags)
  os.chdir(paths.basePath)
  # seperate the the input path from the experiment name
  experimentPath, experimentName = paths.getPathAndName(experimentPathName)
  outscript=experimentName+".run"
  expname = get_EXPNAME(outscript)
  make_runscript_command=paths.basePath+"/config/make_target_runscript "
  inoutFiles="in_script="+experimentPathName+" in_script=exec.iconrun out_script="+outscript+" EXPNAME="+expname+" "
  print(make_runscript_command+inoutFiles+runflags)
  status = runCommand(make_runscript_command+inoutFiles+runflags)
  if not status == 0:
    print("make_runscripts failed")
  os.chdir(paths.thisPath)
  return status, outscript

def get_EXPNAME(runscript):
  expname = os.path.splitext(runscript)[0]
  expname = expname.split(".",1)[1]
  return expname
  

#-----------------------------------------------------------------------

