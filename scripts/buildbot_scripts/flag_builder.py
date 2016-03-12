#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create an experiment list
#==============================================================================
from buildbot_builders import *
from model_paths import *
import argparse

parser = argparse.ArgumentParser(description='List experiments an experiment list.')
parser.add_argument('flags', type=str, help='the flags for the builder')
parser.add_argument('--list', dest="list_name", type=str, help='buildbot list', required=True)
parser.add_argument('--builders', dest="builder_names", nargs='*', help='buildbers', required=True)
args = parser.parse_args()

#print(args.name)
paths = model_paths()
if not paths.thisListExists(args.list_name):
  print("This list does not exist.")
  quit()
  
thisList  = buildbot_experimentList(args.list_name)
thisList.set_builders_flags(args.builder_names, args.flags)
thisList.write()

