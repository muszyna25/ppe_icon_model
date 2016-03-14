#!/usr/bin/python
# -*- coding: utf-8 -*-
#==============================================================================
# driver for testing the builder classes
#==============================================================================
from buildbot_builders import *
from model_paths import *

listname="icon-dev"
paths = model_paths()
thisList  = buildbot_experimentList(listname)
builder_flags, configure_flags, experimentList = thisList.getBuilderProperties("THUNDER_nag")

print(builder_flags)
print(configure_flags)
print(experimentList)

quit()


