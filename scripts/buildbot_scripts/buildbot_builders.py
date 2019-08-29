#! /usr/bin/env python
# -*- coding: utf-8 -*-
#==============================================================================
# Simple library for the builders and experiment lists
#
#
#                            list
#                            /  \
#                           /    \
#             {queue} + machines  experiments
#                         /        /
#                        /        /
#                       /        /
#{configure flags} + builders   /
#{builder flags}        \      /
#                        \    /
#                     associations + {run flags}
#
#
# The class structure is:
#   buildbot_experiments_list
#     buildbot_machines_list
#       buildbot_machines
#         buildbot_builders
#     experimentList
#       experiment

#==============================================================================

import weakref
import sys
import os
from build_interfaces import make_runscript, make_binaries_interface
from model_paths import paths

verbal=True
#paths=None
#-----------------------------------------------------------------------
# contains a list of machines+builders and a list of experiments for them

class buildbot_experiments_list(object):
  
  def __init__(self, name):
    #global paths
    self.name  = name
    self.experimentList = {}
    #paths = model_paths()
    self.buildbot_machines_list = buildbot_machines_list(name)
    if paths.thisListExists(self.name):
      self.read()
    #else:  
      #self.create_all_builders(name) # creates self.buildbot_machines_list with all machines+builders

  # Add an experiment only to the list; no builder is associated here.
  # If the experiment name exists, it just returns the pointer.
  # If one needs to check first the existence of an experiment,
  # should use the get_experiment_value()
  def add_experiment(self, name):
    if not self.experimentList.get(name):
      self.experimentList[name] = buildbot_experiment(name, self)
    return self.experimentList[name]
        
  def add_experimentsByNameToBuilders(self, experimentNames, builder_list, runflags):
    exp_list = []
    for name in experimentNames:
      experiment = self.add_experiment(name)
      exp_list.append(experiment)
      experiment.add_builders(builder_list, runflags, "Normal")
    return exp_list
    
  def add_experimentsByNameToBuildersByName(self, experimentNames, builder_name_list, runflags):    
    builder_list = self.buildbot_machines_list.get_buildersByName(builder_name_list)
    return self.add_experimentsByNameToBuilders(experimentNames, builder_list, runflags)

  def add_experimentsByNameToAllBuilders(self, experimentNames, runflags):
    all_builders = self.buildbot_machines_list.get_all_builders()
    return self.add_experimentsByNameToBuilders( experimentNames, all_builders, runflags)

  def add_experimentsByNameToBuildersWithOptions(self, experimentNames, machinesNames, withFlags, withoutFlags, runflags):
    if machinesNames:
      buildersList = self.buildbot_machines_list.get_machinesBuilders(machinesNames)
    else:
      buildersList = self.buildbot_machines_list.get_all_builders()
    buildersList_withOptions = self.buildbot_machines_list.get_buildersWithOptions(buildersList, withFlags, withoutFlags)
    return self.add_experimentsByNameToBuilders( experimentNames, buildersList_withOptions, runflags)   
    
  def delete_experiment(self, experiment):
    self.experimentList[experiment.name].delete()
    
  # Note: this is only be called by an experiment object; use the delete_experiment() in all other cases
  def delete_experimentFromBuildbotList(self, experiment):
    print("Deleting "+experiment.name+" from "+self.name+" list...")
    del self.experimentList[experiment.name]
    
  def delete_experimentsByName(self, experiment_list):
    for experimentName in experiment_list:
      experiment = self.get_experiment(experimentName)
      self.delete_experiment(experiment)
      
  def delete_experimentsByName_fromBuilders(self, experimentsNames, builders_list):
    for experimentName in experimentsNames:
      experiment = self.get_experiment(experimentName)
      experiment.delete_fromBuilders(builders_list)

  def delete_experimentsByNameFromBuildersByName(self, experiment_list, buildersName):
    builders_list = self.buildbot_machines_list.get_buildersByName(buildersName)
    self.delete_experimentsByName_fromBuilders(experiment_list, builders_list)

  def delete_experimentsByNameFromBuildersWithOptions(self, experimentNames, machinesNames, withFlags, withoutFlags):
    if machinesNames:
      buildersList = self.buildbot_machines_list.get_machinesBuilders(machinesNames)
    else:
      buildersList = self.buildbot_machines_list.get_all_builders()
    buildersList_withOptions = self.buildbot_machines_list.get_buildersWithOptions(buildersList, withFlags, withoutFlags)
    return self.delete_experimentsByName_fromBuilders( experimentNames, buildersList_withOptions)

  def set_builders_flags(self, builders_names, flag):
    builders = self.get_BuildersByName(builders_names)
    for builder in builders:
      builder.set_builder_flags(flag)

  def get_builder_flags(self, builder_name):
    builder = self.get_BuildersByName([builder_name])[0]
    return builder.get_builder_flags()
      
  def set_configure_flags(self, builders_names, flag):
    builders = self.get_BuildersByName(builders_names)
    for builder in builders:
      builder.set_configure_flags(flag)
    
  def get_experiment(self, name):
    experiment = self.experimentList.get(name)
    if not experiment:
      print("Error: experiment "+name+" not found in "+self.name+" list. Stop")
      quit(1)
    return experiment

  def get_experiment_value(self, name):
    experiment = self.experimentList.get(name)
    return experiment
    
  def get_MachineByName(self, name):
    return self.buildbot_machines_list.get_MachineByName(name)
    
  def get_BuildersByName(self, names):
    if names[0] == "*":
      return self.buildbot_machines_list.get_all_builders()
    else:
      return self.buildbot_machines_list.get_buildersByName(names)
        
  def get_BuilderExperimentNames(self, builder_name):
    return self.buildbot_machines_list.get_BuilderExperimentNames(builder_name)
    
  def print_ExperimentsBuilders(self, experimentsNames):
    for experimentName in experimentsNames:
      experiment = self.get_experiment(experimentName)
      experiment.print_builders()

  def print_list(self):
    self.buildbot_machines_list.print_builders()
     
  def create_all_builders(self, name):
    self.buildbot_machines_list.create_all_builders()

  def make_binaries(self, builder_name):
    builder = self.get_BuildersByName([builder_name])[0]
    return builder.make_binaries()

  # if succesful returns a list of the runscripts
  #  otherwise returns the status
  def make_runscripts(self, builder_name):
    builder = self.get_BuildersByName([builder_name])[0]
    return builder.make_runscripts()

  #  returns:
  #    the builder flags (active, build_only, restricted): string
  #    the configure flags: string
  #    the list of expriments to run in a list of the form [[path, name],[path,name],..]
  #         note: the experiment path is the one under the model_paths.runPath
  #
  def get_BuilderProperties(self, builder_name):
    experimentList = []
    builder = self.get_BuildersByName([builder_name])[0]
    #experimentPathNames = builder.getExperimentNames()
    #for experimentPathName in experimentPathNames:
      ## seperate the the input path from the experiment name
      #experimentPath, experimentName = paths.getPathAndName(experimentPathName)
      #experimentList.append([experimentPath, experimentName])
    return builder.get_builder_flags(), builder.get_configure_flags(), builder.getExperimentNames()

  def add_machine(self, name, queue):
    machine=self.buildbot_machines_list.add_machine(name, queue)

  def add_builder(self, builder_name, machine_name, configure_flags, builder_flags):
    machine = self.get_MachineByName(machine_name)
    return machine.add_builder(builder_name,  configure_flags, builder_flags)
    
  #---------------------------------------
  # i/o routines
  #---------------------------------------  
  # the files have a simple ascii form:
  # <keyword>:<name> <parameters>
  # keyword=machine|experiment|builder
  # the hierarchy is:
  # machine:
  #   builder:
  #     experiment:
  def write(self):
    fileName=paths.get_thisListPath(self.name)
    listfile = open(fileName, 'w')
    self.buildbot_machines_list.writeToFile_builders(listfile)
    listfile.close()
    
  def read(self):
    global verbal
    verbal_oldStatus = verbal
    verbal = False
    try:
      fileName=paths.get_thisListPath(self.name)
      listfile = open(fileName, 'r')
      for inLine in listfile:
        inputs  = inLine.rstrip().split('|')
        keyword = inputs[0]
        name    = inputs[1]
        #print(inputs)
        if   (keyword == "machine"):
          machine=self.buildbot_machines_list.add_machine(name, inputs[2])
        elif (keyword == "builder"):
          builder=machine.add_builder(name, inputs[2], inputs[3])
        elif (keyword == "experiment"):
          experiment = self.add_experiment(name)
          # force the addition in case the builder is not active
          experiment.add_builders([builder], inputs[2], "Force")
      listfile.close()     
    except IOError as e:
      print("I/O error({0}): {1}".format(e.errno, e.strerror)+" in reading list "+self.name+". Stop.")
      quit(1)
      
    verbal = verbal_oldStatus
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# a list of machines
class buildbot_machines_list(object):

  def __init__(self, name):
    self.name  = name
    self.machines = {}
    self.builders = weakref.WeakValueDictionary() # {}

  def add_machine(self, name, queue):
    self.machines[name] = buildbot_machine(name, queue, self)
    return self.machines[name]

  # this is only called by a child machine
  def add_builder(self, builder):
    if self.builders.get(builder.name):
      print("Error: trying to add existing builder "+builder.name+". Stop")
      quit(1)
    self.builders[builder.name] = builder
    
  def print_builders(self):
    print("----------------------------")
    sortKeys = sorted(self.machines.keys())
    for machineKey in sortKeys:
      self.machines[machineKey].print_builders()
            
  def writeToFile_builders(self, listfile):
    for machine in self.machines.values():
      machine.writeToFile_builders(listfile)

  def get_machinesBuilders(self, machine_names):
    builders_list = []
    for machine_name in machine_names:
      machine = self.get_MachineByName(machine_name)
      builders_list.extend(machine.get_all_builders())
    return builders_list
    
  def get_all_builders(self):
    return list(self.builders.values())

  def get_buildersWithOptions(self, buildersList, withFlags, withoutFlags):
    buildersWithOptions = []
    for builder in buildersList:
      if builder.hasOptions(withFlags, withoutFlags):
        buildersWithOptions.append(builder)
    return buildersWithOptions

  def get_buildersByName(self,buildersNames):
    buildersList = []
    for builderName in buildersNames:
      builder = self.builders.get(builderName)
      if not builder:
        print("Error: builder "+builderName+" not found in "+self.name+".")
        quit(1)
      else:
        buildersList.append(builder)
    return buildersList
    
  def get_MachineByName(self,machineName):
    machine = self.machines.get(machineName)
    if not machine:
      print("Error: machine "+machineName+" not found in "+self.name+".")
      quit(1)
    return machine

  def get_BuilderExperimentNames(self, builder_name):
    if not self.builders.get(builder_name):
      print("Error: get_BuilderExperimentNames: not existing builder "+builder_name+". Stop")
      quit(1)
    return self.builders[builder_name].getExperimentNames()

  # updates the configure flags based on the builder flags
  #def update_builder_configuration(self):
    #for builder in self.builders.values():
      #builder.update_builder_configuration()
     
  # not used, see create_all_builders script
  #def create_all_builders(self):
    ## add_machine(name, queue)
    ## add_builder(name, configure flags, builder flags)
    ## mistral builders
    #mistral               = self.add_machine('mistral', 'compute')
    #mistral_gcc           = mistral.add_builder('MISTRAL_gcc', '--with-fortran=gcc', 'Active')
    #mistral_intel         = mistral.add_builder('MISTRAL_intel', '--with-fortran=intel', 'Active')
    #mistral_intel_hybrid  = mistral.add_builder('MISTRAL_intel_hybrid', '--with-fortran=intel --with-openmp', 'Active')
    #mistral_intel_openmp  = mistral.add_builder('MISTRAL_intel_openmp', '--with-fortran=intel --without-mpi --with-openmp --without-yac', 'Active')
    #mistral_nag           = mistral.add_builder('MISTRAL_nag', '--with-fortran=nag', 'Active')
    #mistral_ocean         = mistral.add_builder('MISTRAL_ocean', '--with-fortran=intel --with-openmp', 'Ocean')
    #mistral_ocean         = mistral.add_builder('MISTRAL_ocean', '--with-fortran=intel --with-openmp --with-flags=ocean', 'Ocean')
## CSCS builders
    #daint                 = self.add_machine('daint', 'default')
    #daint_cpu_cce         = daint.add_builder('DAINT_CPU_cce', '', 'Active')
## breeze builders
    #breeze                = self.add_machine('breeze', 'default')
    #breeze_gcc            = breeze.add_builder('BREEZE_gcc', '--with-fortran=gcc', 'build_only')
    #breeze_gcc_openmp     = breeze.add_builder('BREEZE_gcc_openmp', '--with-fortran=gcc --with-openmp', 'build_only')
    #breeze_intel          = breeze.add_builder('BREEZE_intel', '--with-fortran=intel', 'build_only')
    #breeze_intel_openmp   = breeze.add_builder('BREEZE_intel_openmp', '--with-fortran=intel --with-openmp', 'build_only')
    #breeze_nag            = breeze.add_builder('BREEZE_nag', '--with-fortran=nag', 'build_only')
## DWD builders
    #dwd                   = self.add_machine('dwd', 'default')
    #dwd_cray              = dwd.add_builder('DWD_cray', '', 'Active')
    #dwd_cray_production   = dwd.add_builder('DWD_cray_production', '', 'Active')

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# one machine; holds a number of builders, as well as the queue to be
# used for the associated list
class buildbot_machine(object):

  def __init__(self, name, queue, machine_list):
    self.name  = name
    self.queue = queue
    self.machine_list = machine_list
    self.builders= {}

  def add_builder(self, name, configure_flags, builder_flags):
    self.builders[name] = buildbot_builder(name, self, configure_flags, builder_flags)
    self.machine_list.add_builder(self.builders[name])
    return self.builders[name]

  def print_builders(self):
    print(self.name+" ("+self.queue+"):")
    sortKeys = sorted(self.builders.keys())
    for builderKey in sortKeys:
      self.builders[builderKey].print_builder_experiments()

  def get_all_builders(self):
    return list(self.builders.values())

  def writeToFile_builders(self, listfile):
    listfile.write("machine|"+self.name+"|"+self.queue+"\n")
    for builder in self.builders.values():
      builder.writeToFile_builder_experiments(listfile)

#-----------------------------------------------------------------------

      
#-----------------------------------------------------------------------
# experiment class
class buildbot_experiment(object):

  def __init__(self, name, experimentList):
    self.name  = name
    self.builders = weakref.WeakValueDictionary() # {}
    self.builder_runflags = {}
    self.experimentList = experimentList # weakref.ref(experimentList)

  def add_builders(self, buildersList, runflags, ActionFlag):
    for builder in buildersList:
      if builder.isActive() or ActionFlag == "Force":
        if verbal: print("adding "+self.name+" to "+builder.name+"...")
        self.builders[builder.name] = builder
        self.builder_runflags[builder.name] = runflags
        builder.add_experiment_onlyFromExperimentObject(self, runflags)
    
  def print_experiment(self):
    print("    "+self.name)
    
  def writeToFile_experiment(self,listfile, builder_name):
    listfile.write("experiment|"+self.name+"|"+self.builder_runflags[builder_name]+"\n")
    
  def print_builders(self):
    print("----------------------------")
    print("Experiment:"+self.name+"; builders:")
    sortKeys = sorted(self.builders.keys())
    for builderKey in sortKeys:
      print("   "+builderKey)   
    #for builder in self.builders.values():
      #builder.print_builder()

  def delete(self):
    print("Deleting experiment "+self.name+" from all builds...")
    for builder in self.builders.values():
      builder.delete_experiment_onlyFromExperimentObject(self)
    self.builders.clear()
    self.experimentList.delete_experimentFromBuildbotList(self)
    
  def delete_fromBuilders(self, builders):
    for builder in builders:
      if self.builders.get(builder.name):
        print("Deleting experiment "+self.name+" from builder "+builder.name+"...")
        builder.delete_experiment_onlyFromExperimentObject(self)
        del self.builders[builder.name]
      
  #def __del__(self):
    #print("Deleting experiment "+self.name+" done.")
#-----------------------------------------------------------------------
   

#-----------------------------------------------------------------------
# the builder class; belongs to a machine and stores the flags associated
# with it. It also holds a list of experiments associated with the builder;
# this list is driven through the experiment class, not directly from this class
class buildbot_builder(object):

  def __init__(self, name, machine, configure_flags, builder_flags):
    self.name  = name
    self.machine = machine
    self.configure_flags = configure_flags
    self.builder_flags = builder_flags
    self.experiments = weakref.WeakValueDictionary() # {}
    self.experiments_runflags = {} # {}

  def isActive(self):
    return self.builder_flags == "Active"
    
  def set_builder_flags(self, flags):
    self.builder_flags = flags

  def get_builder_flags(self):
    return self.builder_flags

  def set_configure_flags(self, flags):
    self.configure_flags = flags

  def get_configure_flags(self):
    return self.configure_flags
    
  # updates the configure flags based on the builder flags
  #def update_builder_configuration(self):
    #ocean_flags=["--disable-atmo","--disable-jsbach","--with-yac=no","--with-flags=ocean"]
    #aes_flags=["--disable-ocean"]
    #if "Ocean" in self.builder_flags:
      ## add  ocean_flags to configure flags
      #for flag in ocean_flags:
        #if not flag in self.configure_flags:
          #self.configure_flags+=" "+flag
         
  def hasOptions(self, withFlags, withoutFlags):
    hasTheseOptions = True
    if withFlags:
      for flag in withFlags:
        hasTheseOptions = hasTheseOptions and flag in self.configure_flags
    if withoutFlags:
      for flag in withoutFlags:
        hasTheseOptions = hasTheseOptions and not flag in self.configure_flags      
    return hasTheseOptions
    
  def getExperimentNames(self):
    return list(self.experiments.keys())
  
  def getExperiments(self):
    return list(self.experiments.values())

  def print_builder(self):
    print("  "+self.name+" ("+self.builder_flags+"):"+self.configure_flags)

  # this should be called only from an experiment object
  def add_experiment_onlyFromExperimentObject(self, experiment, runflags):
    self.experiments[experiment.name] = experiment
    self.experiments_runflags[experiment.name] = runflags
    
  # this should be called only from an experiment object
  def delete_experiment_onlyFromExperimentObject(self, experiment):
    print(" Deleting "+experiment.name+" from "+self.name+"...")
    del self.experiments[experiment.name]

  def make_binaries(self):
    return make_binaries_interface(self.configure_flags, self.builder_flags)

  # if succesful returns a list of the runscripts
  #  otherwise returns the status
  def make_runscripts(self):
    os.chdir(paths.basePath)
    runscriptList = []
    if "build_only" in self.builder_flags or "Inactive" in self.builder_flags :
      return runscriptList
    for experiment in self.experiments.values():
      experimentPathName = experiment.name
      status, runscript = make_runscript(experimentPathName, self.experiments_runflags[experimentPathName])
      if not status == 0:
        return status
      runscriptList.append(runscript)
    return runscriptList

  def print_builder_experiments(self):
    self.print_builder()
    sortKeys = sorted(self.experiments.keys())
    for experimentKey in sortKeys:
      print("    "+experimentKey+" runflags:"+self.experiments_runflags[experimentKey])
      
    #for experiment in self.experiments.values():
      #experiment.print_experiment()

  def writeToFile_builder_experiments(self, listfile):
    listfile.write("builder|"+self.name+'|'+self.configure_flags+"|"+self.builder_flags+"\n")
    for experiment in self.experiments.values():
      experiment.writeToFile_experiment(listfile, self.name)
#-----------------------------------------------------------------------

