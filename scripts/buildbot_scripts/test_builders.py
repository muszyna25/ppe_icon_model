#! /usr/bin/env python
# -*- coding: utf-8 -*-
#==============================================================================
# driver for testing the builder classes
#==============================================================================
from buildbot_builders import *
from model_paths import *

listname="icon-dev"
paths = model_paths()
thisList  = buildbot_experiments_list(listname)
builder_flags, configure_flags, experimentList = thisList.getBuilderProperties("THUNDER_gcc")

print(builder_flags)
print(configure_flags)
print(experimentList)

#thisList.make_binaries("THUNDER_gcc")
runscriptsList = thisList.make_runscripts("THUNDER_gcc")
print(runscriptsList)

quit()


