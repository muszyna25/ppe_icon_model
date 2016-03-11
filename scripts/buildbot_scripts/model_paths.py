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
    self.experimentsListPath = self.thisPath+"/experiment_lists"
    
  def get_experimentsNames_inPath(self, pathInRun):
    return glob.glob(self.runPath+"/"+pathInRun)

  def get_thisExperimentListPath(self, listName):
    return self.experimentsListPath+"/"+listName
    
  def thisExperimentListExists(self, listName):
    return os.path.isfile(self.get_thisExperimentListPath(listName))
    
  def deleteThisExperimentList(self, listName):
    if not self.thisExperimentListExists(listName):
      print("The list "+listName+" does not exist.")
      quit()
    os.remove(self.get_thisExperimentListPath(listName))


  def print_paths(self):
    print("Base path:"+self.basePath)
    print("Run path:"+self.runPath)
    print("This path:"+self.thisPath)
    

