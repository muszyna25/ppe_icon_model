#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create an experiment list
#==============================================================================
from buildbot_builders import *
from model_paths import *
import argparse

parser = argparse.ArgumentParser(description='List experiments an experiment list.')
parser.add_argument('name', type=str, help='the name of the list')
args = parser.parse_args()

#print(args.name)
name=args.name
paths = model_paths()
if not paths.thisListExists(name):
  print("This list does not exist.")
  quit()
  
thisList  = buildbot_experimentList(name)
thisList.print_list()

