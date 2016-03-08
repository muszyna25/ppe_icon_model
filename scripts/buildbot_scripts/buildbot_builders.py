#!/usr/bin/python
# -*- coding: utf-8 -*-
#==============================================================================
# Simple lib for the builders
#==============================================================================




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
    for name, builder in self.builders.iteritems():
      #print name, builder.name
      builder.print_builder()


#-----------------------------------------------------------------------
class buildbot_builder(object):

  def __init__(self, name, machine, flags):
    self.name  = name
    self.machine = machine
    self.configure_flags = flags
    self.experiments = {}

  def print_builder(self):
    print self.name+':', self.configure_flags

  def add_experiment(self, experiment):
    self.experiments[name] = experiment
    
#-----------------------------------------------------------------------
class buildbot_experiment(object):
  
  def __init__(self, name):
    self.name  = name
    self.builders = {}

  def add_builder(self, builder):
    self.builders[builder.name] = builder
    
  
#-----------------------------------------------------------------------
class buildbot_machine_list(object):
  
  def __init__(self, name):
    self.name  = name
    self.machines = {}

  def add(self, name, queue):
    self.machines[name] = buildbot_machine(name, queue)
    return self.machines[name]
    
  def print_builders(self):
    for name, machine in self.machines.iteritems():
      print name+":"
      machine.print_builders()
#-----------------------------------------------------------------------

# main
def create_all_builders(name):
  buildobot_machines    =  buildbot_machine_list(name)
  
  mistral               = buildobot_machines.add('mistral', 'compute')
  mistral_gcc           = mistral.add_builder('mistral_gcc', '--with-fortran=gcc')
  mistral_intel         = mistral.add_builder('mistral_intel', '--with-fortran=intel')
  mistral_intel_hybrid  = mistral.add_builder('mistral_intel_hybrid', '--with-fortran=intel --with-openmp')
  mistral_intel_openmp  = mistral.add_builder('mistral_intel_openmp', '--with-fortran=intel --without-mpi --with-openmp --without-yac ')
  mistral_nag           = mistral.add_builder('mistral_nag', '--with-fortran=nag')
  mistral_nag_mtime     = mistral.add_builder('mistral_nag_mtime', '--with-fortran=nag')
  mistral_nag_serial    = mistral.add_builder('mistral_nag_serial', '--with-fortran=nag --without-mpi')

  return buildobot_machines
#-----------------------------------------------------------------------

  
buildobot_machines = create_all_builders("buildbot")

buildobot_machines.print_builders()


