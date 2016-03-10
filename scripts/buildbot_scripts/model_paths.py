#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# model paths classes
#==============================================================================
import os
import glob

#-----------------------------------------------------------------------
class model_paths(object):
  
  def __init__(self):
    self.thisPath   = os.getcwd()
    splitPath       = os.path.split(self.thisPath)
    splitPath       = os.path.split(splitPath[0])
    self.basePath   = splitPath[0]
    self.runPath    = self.basePath+"/run"
    
  def get_experimentsNames_inPath(self, pathInRun):
    return glob.glob(self.runPath+"/"+pathInRun)

  def print_paths(self):
    print("Base path:"+self.basePath)
    print("Run path:"+self.runPath)
    print("This path:"+self.thisPath)
    

