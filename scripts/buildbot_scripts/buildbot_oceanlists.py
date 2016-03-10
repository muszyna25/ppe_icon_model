#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create the ocean experiment lists
#==============================================================================
from buildbot_builders import *
import os

thisPath = os.getcwd()
print(thisPath)

#quit()


oceanList  = buildbot_experimentList("ocean")

myexp      = oceanList.add_experimentsByNameToAllBuildersWithoutFlag(["checksuite.ocean_internal/omip/exp.ocean_omip_ptest"], "without-mpi")
myexp      = oceanList.add_experimentsByNameToAllBuilders(["checksuite.ocean_internal/omip/exp.ocean_omip_short"])
myexp      = oceanList.add_experimentsByNameToAllBuildersWithoutFlag(["checksuite.ocean_internal/technical/exp.test_ocean_omip_technical"], "without-mpi")

oceanList.print_list()

oceanList.write()

sameOceanList  = buildbot_experimentList("ocean")

sameOceanList.read()
sameOceanList.print_list()

myexp = sameOceanList.get_experiment("checksuite.ocean_internal/omip/exp.ocean_omip_short")
myexp.print_builders()
sameOceanList.delete_experiment(myexp)

# this will cause an error, since this experiment is deleted
myexp = sameOceanList.get_experiment("checksuite.ocean_internal/omip/exp.ocean_omip_short")


