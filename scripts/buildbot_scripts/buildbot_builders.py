#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# Simple library for the builders and experiment lists
#
# The class structure is:
#   buildbot_experimentList
#     buildbot_machine_list
#       buildbot_machines
#         buildbot_builders
#     experimentList
#       experiment
#==============================================================================

import weakref
import sys
from model_paths import *

#-----------------------------------------------------------------------
# contains a list of machines+builders and a list of experiments for them
class buildbot_experimentList(object):

  def __init__(self, name):
    self.name  = name
    self.buildbot_machine_list = self.create_all_builders(name)
    self.experimentList = {}
    self.paths = model_paths()

  # Add an experiment only to the list; no builder is associated here.
  # If the experiment name exists, it just returns the pointer.
  # If one needs to check first the existence of an experiment,
  # should use the get_experiment_state()
  def add_experiment(self, name):
    if not self.experimentList.get(name):
      self.experimentList[name] = buildbot_experiment(name, self)
    return self.experimentList[name]

  def add_experimentByNameToBuilder(self, name, builder):
    experiment = self.add_experiment(name)
    experiment.add_builder(builder)
    return experiment
    
  def add_experimentsByNameToMachine(self, name_list, machine_name):
    builders_list = self.buildbot_machine_list.get_machine_builders(machine_name)
    self.add_experimentsByNameToBuilders(name_list, builders_list)
    
  def add_experimentsByNameToBuilders(self, name_list, builder_list):
    exp_list = []
    for name in name_list:
      for builder in builder_list:
        exp_list.append(self.add_experimentByNameToBuilder(name, builder))
    return exp_list
    
  def add_experimentsByNameToBuildersByName(self, name_list, builder_name_list):
    builder_list = []
    for builder_name in builder_name_list:
      builder_list.append(sel.getBuilderByName(builder_name))
    return add_experimentsByNameToBuilders(self, name_list, builder_list)

  def add_experimentsByNameToAllBuilders(self, name_list):
    all_builders = self.buildbot_machine_list.get_all_builders()
    return self.add_experimentsByNameToBuilders( name_list, all_builders)
    
  def add_experimentsByNameToBuildersWithFlag(self, name_list, flag):
    builders_list = self.buildbot_machine_list.get_builders_withFlag(flag)
    return self.add_experimentsByNameToBuilders( name_list, builders_list)
    
  def add_experimentsByNameToBuildersWithoutFlag(self, name_list, flag):
    builders_list = self.buildbot_machine_list.get_builders_withoutFlag(flag)
    return self.add_experimentsByNameToBuilders( name_list, builders_list)
    
  def delete_experiment(self, experiment):
    self.experimentList[experiment.name].delete()

  # Note: this should only be called by an experiment object; use the delete_experiment() in all other cases
  def delete_experimentFromBuildbotList(self, experiment):
    print("Deleting "+experiment.name+" from "+self.name+" list...")
    del self.experimentList[experiment.name]
    
  def delete_experimentsByName(self, experiment_list):
    for experimentName in experiment_list:
      experiment = self.get_experiment(experimentName)
      self.delete_experiment(experiment)
      
  def delete_experimentsByName_fromBuilders(self, experiment_list, builders_list):
    for experimentName in experiment_list:
      experiment = self.get_experiment(experimentName)
      for builder in builders_list:
        experiment.delete_fromBuilder(builder)

  def delete_experimentsByName_fromBuildersName(self, experiment_list, buildersName):
    builders_list = self.getBuildersListByName(buildersName)
    self.delete_experimentsByName_fromBuilders(experiment_list, builders_list)
            
  def delete_experimentsByName_fromMachineName(self, experiment_list, machine_name):
    builders_list = self.buildbot_machine_list.get_machine_builders(machine_name)
    self.delete_experimentsByName_fromBuilders(experiment_list, builders_list)

  def delete_experimentsByName_fromBuildersWithFlag(self, experiment_list, flag):
    builders_list = self.buildbot_machine_list.get_builders_withFlag(flag)
    self.delete_experimentsByName_fromBuilders(experiment_list, builders_list)

  def delete_experimentsByName_fromBuildersWithoutFlag(self, experiment_list, flag):
    builders_list = self.buildbot_machine_list.get_builders_withoutFlag(flag)
    self.delete_experimentsByName_fromBuilders(experiment_list, builders_list)
        
  def get_experiment(self, name):
    experiment = self.experimentList.get(name)
    if not experiment:
      print("Error: experiment "+name+" not found in "+self.name+" list. Stop")
      quit()
    return experiment

  def get_experiment_state(self, name):
    experiment = self.experimentList.get(name)
    return experiment
    
  def getMachineByName(self, name):
    return self.buildbot_machine_list.getMachineByName(name)

  def getBuildersListByName(self, name_list):
    builders_list = []
    for name in name_list:
      builders_list.append(self.getBuilderByName(name))
    return builders_list
    
  def getBuilderByName(self, name):
    return self.buildbot_machine_list.getBuilderByName(name)
        
  def getBuilderExperimentNames(self, builder_name):
    return self.buildbot_machine_list.getBuilderExperimentNames(builder_name)

  def print_list(self):
    self.buildbot_machine_list.print_builders()
     
  def create_all_builders(self, name):
    buildbot_machines     =  buildbot_machine_list(name)
    # mistral builders
    mistral               = buildbot_machines.add_machine('mistral', 'compute')
    mistral_gcc           = mistral.add_builder('mistral_gcc', '--with-fortran=gcc')
    mistral_intel         = mistral.add_builder('mistral_intel', '--with-fortran=intel')
    mistral_intel_hybrid  = mistral.add_builder('mistral_intel_hybrid', '--with-fortran=intel --with-openmp')
    mistral_intel_openmp  = mistral.add_builder('mistral_intel_openmp', '--with-fortran=intel --without-mpi --with-openmp --without-yac ')
    mistral_nag           = mistral.add_builder('mistral_nag', '--with-fortran=nag')
    mistral_nag_mtime     = mistral.add_builder('mistral_nag_mtime', '--with-fortran=nag')
    mistral_nag_serial    = mistral.add_builder('mistral_nag_serial', '--with-fortran=nag --without-mpi')
    # CSCS builders
    daint_cpu             = buildbot_machines.add_machine('daint_cpu', 'default')
    daint_cpu_cce         = daint_cpu.add_builder('DAINT_CPU_cce', '')
    # thunder builders
    thunder               = buildbot_machines.add_machine('thunder', 'mpi-compute')
    thunder_gcc           = thunder.add_builder('thunder_gcc', '--with-fortran=gcc')
    thunder_intel_hybrid  = thunder.add_builder('thunder_intel_hybrid', '--with-fortran=intel --with-openmp')
    thunder_nag           = thunder.add_builder('thunder_nag', '--with-fortran=nag')
         
    return buildbot_machines

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
    fileName=self.paths.get_thisExperimentListPath(self.name)
    listfile = open(fileName, 'w')
    self.buildbot_machine_list.writeToFile_builders(listfile)
    listfile.close()
    
  def read(self):
    try:
      fileName=self.paths.get_thisExperimentListPath(self.name)
      listfile = open(fileName, 'r')
      for inLine in listfile:
        inputs=inLine.split(':')
        keyword=inputs[0]
        parameters=inputs[1].rstrip().split(' ',2)
        name=parameters[0]
        #print(keyword, name)
        if   (keyword == "machine"):
          machine=self.getMachineByName(name)
        elif (keyword == "builder"):
          builder=self.getBuilderByName(name)
        elif (keyword == "experiment"):
          self.add_experimentByNameToBuilder(name,builder)

      listfile.close()
      
    except IOError as e:
      print("I/O error({0}): {1}".format(e.errno, e.strerror)+" in reading list "+self.name+". Stop.")
      quit()
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# a list of machines
class buildbot_machine_list(object):

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
      quit()
    self.builders[builder.name] = builder
    
  def print_builders(self):
    print("----------------------------")
    for machine in self.machines.values():
      machine.print_builders()
      
  def writeToFile_builders(self, listfile):
    for machine in self.machines.values():
      machine.writeToFile_builders(listfile)

  def get_machine_builders(self, machine_name):
    return self.machines[machine_name].get_all_builders()
    
  def get_all_builders(self):
    all_builders = []
    for machine in self.machines.values():
      all_builders.extend(machine.get_all_builders())
    return all_builders
    
  def get_builders_withFlag(self, flag):
    builders = []
    for machine in self.machines.values():
      builders.extend(machine.get_builders_withFlag(flag))
    return builders

  def get_builders_withoutFlag(self, flag):
    builders = []
    for machine in self.machines.values():
      builders.extend(machine.get_builders_withoutFlag(flag))
    return builders

  def getBuilderByName(self,builderName):
    for machine in self.machines.values():
      builder = machine.getBuilderByName(builderName)
      if builder:
        return builder
    print("Error: builder "+builderName+" not found in "+self.name+".")
    quit()
    
  def getMachineByName(self,machineName):
    machine = self.machines.get(machineName)
    if not machine:
      print("Error: machine "+machineName+" not found in "+self.name+".")
      quit()

  def getBuilderExperimentNames(self, builder_name):
    if not self.builders.get(builder_name):
      print("Error: getBuilderExperimentNames: not existing builder "+builder_name+". Stop")
      quit()
    return self.builders[builder_name].getExperimentNames()
     
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

  def add_builder(self, name, flags):
    self.builders[name] = buildbot_builder(name, self, flags)
    self.machine_list.add_builder(self.builders[name])
    return self.builders[name]

  def print_builders(self):
    print(self.name+":")
    for builder in self.builders.values():
      builder.print_builder_experiments()

  def get_all_builders(self):
    return self.builders.values()

  def get_builders_withFlag(self, flag):
    builders = []
    for builder in self.builders.values():
      if builder.hasFlag(flag):
        builders.append(builder)
    return builders

  def get_builders_withoutFlag(self, flag):
    builders = []
    for builder in self.builders.values():
      if not builder.hasFlag(flag):
        builders.append(builder)
    return builders

  def writeToFile_builders(self, listfile):
    listfile.write("machine:"+self.name+" "+self.queue+"\n")
    for builder in self.builders.values():
      builder.writeToFile_builder_experiments(listfile)

  def getBuilderByName(self,name):
    return self.builders.get(name)
#-----------------------------------------------------------------------

      
#-----------------------------------------------------------------------
# experiment class
class buildbot_experiment(object):

  def __init__(self, name, experimentList):
    self.name  = name
    self.builders = weakref.WeakValueDictionary() # {}
    self.experimentList = experimentList # weakref.ref(experimentList)

  def add_builder(self, builder):
    self.builders[builder.name] = builder
    builder.add_experiment_onlyFromExperimentObject(self)

  def add_builderByName(self, builder_name):
    builder = self.experimentList.getBuilderByName(builder_name)
    self.add_builder(builder)
    
  def print_experiment(self):
    print("    "+self.name)
    
  def writeToFile_experiment(self,listfile):
    listfile.write("experiment:"+self.name+"\n")
    
  def print_builders(self):
    print("----------------------------")
    print("Experiment:"+self.name+"; builders:")
    for builder in self.builders.values():
      builder.print_builder()

  def delete(self):
    print("Deleting experiment "+self.name+" from all builds...")
    for builder in self.builders.values():
      builder.delete_experiment_onlyFromExperimentObject(self)
    #print("Deleting experiment "+self.name+" builders...")
    self.builders.clear()
    #print("Deleting experiment "+self.name+" from "+self.experimentList.name+" list...")
    self.experimentList.delete_experimentFromBuildbotList(self)
    #del self
    
  def delete_fromBuilder(self, builder):
    if self.builders.get(builder.name):
      print("Deleting experiment "+self.name+" from builder "+builder.name+"...")
      builder.delete_experiment_onlyFromExperimentObject(self)
      del self.builders[builder.name]
    else:
      print("Warning: delete failed. "+self.name+" experiment is not in the "+ builder.name+" builder.")
      
  #def __del__(self):
    #print("Deleting experiment "+self.name+" done.")
#-----------------------------------------------------------------------
   

#-----------------------------------------------------------------------
# the builder class; belongs to a machine and stores the flags associated
# with it. It also holds a list of experiments associated with the builder;
# this list is driven through the experiment class, not directly from this class
class buildbot_builder(object):

  def __init__(self, name, machine, flags):
    self.name  = name
    self.machine = machine
    self.configure_flags = flags
    self.experiments = weakref.WeakValueDictionary() # {}

  def print_builder(self):
    print("  "+self.name+':', self.configure_flags)

  def hasFlag(self,flag):
    return (flag in self.configure_flags)
  
  def print_builder_experiments(self):
    print("  "+self.name+':', self.configure_flags)
    for experiment in self.experiments.values():
      experiment.print_experiment()

  def getExperimentNames(self):
    return list(self.experiments.keys())
      
  def writeToFile_builder_experiments(self, listfile):
    listfile.write("builder:"+self.name+" "+self.configure_flags+"\n")
    for experiment in self.experiments.values():
      experiment.writeToFile_experiment(listfile)

  # this should be called only from an experiment object
  def add_experiment_onlyFromExperimentObject(self, experiment):
    self.experiments[experiment.name] = experiment
    
  # this should be called only from an experiment object
  def delete_experiment_onlyFromExperimentObject(self, experiment):
    print(" Deleting "+experiment.name+" from "+self.name+"...")
    del self.experiments[experiment.name]
    #print(experiment.name+" is deleted from "+self.name+"...")
#-----------------------------------------------------------------------

