#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create an experiment list
#==============================================================================
from buildbot_builders import *
from model_paths import *
import argparse

parser = argparse.ArgumentParser(description='Create an experiment list.')
parser.add_argument('name', type=str, help='the name of the list')
args = parser.parse_args()

#print(args.name)
name=args.name
paths = model_paths()
if paths.thisExperimentListExists(name):
  print("This list exists.")
  quit()
  
myList  = buildbot_experimentList(name)
myList.write()
print("Experiment list "+name+" is created.")

