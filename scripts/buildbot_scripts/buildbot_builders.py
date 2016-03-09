#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# Simple lib for the builders
#==============================================================================
import weakref

listsFolder="experiment_lists/"
#-----------------------------------------------------------------------
class buildbot_experimentList(object):

  def __init__(self, name):
    self.name  = name
    self.buildbot_machine_list = self.create_all_builders(name)
    self.experimentList = {}

  def add_experiment(self, name):
    if not self.experimentList.get(name):
      self.experimentList[name] = buildbot_experiment(name, self)
    return self.experimentList[name]

  def add_experimentToBuilder(self, name, builder):
    experiment = self.add_experiment(name)
    experiment.add_builder(builder)

  def add_experimentToAllBuilders(self, experiment):
    self.buildbot_machine_list.add_experimentToAllMachines(experiment)

  def add_experimentByNameToAllBuilders(self, name):
    experiment = self.add_experiment(name)
    self.add_experimentToAllBuilders(experiment)
    return experiment

  def delete_experiment(self, experiment):
    self.experimentList[experiment.name].delete()

  # this should only called by an experiment object; use the delete_experiment() in all other cases
  def delete_experimentFromList(self, experiment):
    print("Deleting "+experiment.name+" from "+self.name+" list...")
    del self.experimentList[experiment.name]
    
  def delete_experimentByName(self, experimentName):
    experiment = self.get_experiment(experimentName)
    self.delete_experiment(experiment)

  def delete_experimentByName_fromBuilder(self, experimentName, builderName):
    experiment = self.get_experiment(experimentName)
    experiment.delete_fromBuilder(builderName)
    
  def get_experiment(self, name):
    experiment = self.experimentList.get(name)
    if not experiment:
      print("Error: experiment "+name+" not found in "+self.name+" list. Stop")
      quit()
    return experiment
    
  def getMachineByName(self, name):
    return self.buildbot_machine_list.getMachineByName(name)
    
  def getBuilderByName(self, name):
    return self.buildbot_machine_list.getBuilderByName(name)
        
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
  # <keyword>: name parameters
  # keyword=machine,experiment,builder
  # the hierarchy is:
  # machine:
  #   builder:
  #     experiment:
  def write(self):
    fileName=listsFolder+self.name
    listfile = open(fileName, 'w')
    self.buildbot_machine_list.write_builders(listfile)    
    listfile.close()
    
  def read(self):
    fileName=listsFolder+self.name
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
        self.add_experimentToBuilder(name,builder)
      
    listfile.close()
    
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
class buildbot_machine_list(object):

  def __init__(self, name):
    self.name  = name
    self.machines = {}

  def add_machine(self, name, queue):
    self.machines[name] = buildbot_machine(name, queue)
    return self.machines[name]

  def print_builders(self):
    print("----------------------------")
    for machine in self.machines.values():
      machine.print_builders()
      
  def write_builders(self, listfile):    
    for machine in self.machines.values():
      machine.write_builders(listfile)
      
  def add_experimentToAllMachines(self,experiment):
    for machine in self.machines.values():
      machine.add_experimentToAllBuilders(experiment)

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

      
#-----------------------------------------------------------------------
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
    
  def write_experiment(self,listfile):
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
    print("Deleting experiment "+self.name+" builders...")
    self.builders.clear()
    #print("Deleting experiment "+self.name+" from "+self.experimentList.name+" list...")
    self.experimentList.delete_experimentFromList(self)
    #del self
    
  def delete_fromBuilder(self, builderName):
    print("Deleting experiment "+self.name+" from builder "+builderName+"...")
    self.builders[builderName].delete_experiment_onlyFromExperimentObject(self)
    del self.builders[builderName]
    
  def __del__(self):
    print("Deleting experiment "+self.name+" done.")
   

#-----------------------------------------------------------------------
class buildbot_machine(object):

  def __init__(self, name, queue):
    self.name  = name
    self.queue = queue
    self.builders= {}

  def add_builder(self, name, flags):
    self.builders[name] = buildbot_builder(name, self, flags)
    return self.builders[name]
    
  def print_builders(self):
    print(self.name+":")
    for builder in self.builders.values():
      builder.print_builder_experiments()
      
  def add_experimentToAllBuilders(self,experiment):
    for builder in self.builders.values():
      experiment.add_builder(builder)
      
  def write_builders(self, listfile):
    listfile.write("machine:"+self.name+" "+self.queue+"\n")
    for builder in self.builders.values():
      builder.write_builder_experiments(listfile)

  def getBuilderByName(self,name):
    return self.builders.get(name)


#-----------------------------------------------------------------------
class buildbot_builder(object):

  def __init__(self, name, machine, flags):
    self.name  = name
    self.machine = machine
    self.configure_flags = flags
    self.experiments = weakref.WeakValueDictionary() # {}

  def print_builder(self):
    print("  "+self.name+':', self.configure_flags)
      
  def print_builder_experiments(self):
    print("  "+self.name+':', self.configure_flags)
    for experiment in self.experiments.values():
      experiment.print_experiment()
      
  def write_builder_experiments(self, listfile):
    listfile.write("builder:"+self.name+" "+self.configure_flags+"\n")
    for experiment in self.experiments.values():
      experiment.write_experiment(listfile)

  # this should be called only from an experiment object
  def add_experiment_onlyFromExperimentObject(self, experiment):
    self.experiments[experiment.name] = experiment
    
  # this should be called only from an experiment object
  def delete_experiment_onlyFromExperimentObject(self, experiment):
    print(" Deleting "+experiment.name+" from "+self.name+"...")
    del self.experiments[experiment.name]
    #print(experiment.name+" is deleted from "+self.name+"...")
    
  
#-----------------------------------------------------------------------

