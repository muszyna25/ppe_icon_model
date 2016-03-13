#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create an experiment list
#==============================================================================
from buildbot_builders import *
from model_paths import *
import argparse

parser = argparse.ArgumentParser(description='Add experiments to an experiment list.')
parser.add_argument('experiment_names', type=str, nargs='+', help='experiment names')
parser.add_argument('--list', dest="list_name", type=str, help='buildbot list', required=True)
parser.add_argument('--builders', dest="builder_names", nargs='*', help='buildbers', required=False)
parser.add_argument('--machines', dest="machine_names", nargs='*', help='machine name', required=False)
parser.add_argument('--with-flags', dest="withFlags", nargs='*', help='with build flag', required=False)
parser.add_argument('--without-flags', dest="withoutFlags", nargs='*', help='without build flag', required=False)
parser.add_argument('--runflags', dest="runFlags", type=str, help='use these run flags', required=False)

args = parser.parse_args()

#print(args.experiment_names)
#print(args.list_name)
#print(args.builder_names)
#print(args.machine_names)
#print(args.withFlags)
#print(args.withoutFlags)


if args.builder_names and (args.machine_names or args.withFlags or args.withoutFlags):
  print("You cannot specify a builder with additional options.")
  quit()

paths = model_paths()
experiment_names = paths.get_experimentsNames_inPaths(args.experiment_names)
#print(experiment_names)
if not paths.thisListExists(args.list_name):
  print("The list "+args.list_name+" does not exist.")
  quit()

if not args.runFlags: args.runFlags = ""

thisList  = buildbot_experimentList(args.list_name)

if args.builder_names:
  thisList.add_experimentsByNameToBuildersByName(experiment_names, args.builder_names, args.runFlags)
else:
  thisList.add_experimentsByNameToBuildersWithOptions(experiment_names, args.machine_names, args.withFlags, args.withoutFlags, args.runFlags)

thisList.write()
#thisList.print_list()
print("Experiment list "+args.list_name+" is updated.")

quit()


