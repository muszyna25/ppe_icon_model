#!/usr/bin/python
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
    
  def get_experimentsNames_inPaths(self, pathsInRun):
    os.chdir(self.runPath)
    experimentsNames = []
    for path in pathsInRun:
      experiments = glob.glob(path)
      if not experiments:
        print("No experiment "+path+" found. Stop")
        quit()
      experimentsNames.extend(experiments)
    os.chdir(self.thisPath)
    return experimentsNames

  def get_thisListPath(self, listName):
    return self.experimentsListPath+"/"+listName
    
  def thisListExists(self, listName):
    return os.path.isfile(self.get_thisListPath(listName))
    
  def deleteThisList(self, listName):
    if not self.thisListExists(listName):
      print("The list "+listName+" does not exist.")
      quit()
    os.remove(self.get_thisListPath(listName))

  def print_paths(self):
    print("Base path:"+self.basePath)
    print("Run path:"+self.runPath)
    print("This path:"+self.thisPath)

  def getPathAndName(self, PathName):
    dirName  = os.path.dirname(PathName)
    fileName = os.path.basename(PathName)
    return dirName, fileName
    

