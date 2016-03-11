#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create the ocean experiment lists
#==============================================================================
from buildbot_builders import *
from model_paths import *
  

myPaths = model_paths()
omipFileList = myPaths.get_experimentsNames_inPath("checksuite.ocean_internal/omip/exp.*")
#print(omipFileList)

noList  = buildbot_experimentList("tt")
noList.read()

noList.print_list()
quit()

omipList  = buildbot_experimentList("omip")
myexps    = omipList.add_experimentsByNameToBuildersWithoutFlag(omipFileList, "without-mpi")

omipList.print_list()
omipList.delete_experimentsByName_fromMachineName(omipFileList, "daint_cpu")
omipList.print_list()
omipList.delete_experimentsByName_fromMachineName(omipFileList, "thunder")
omipList.print_list()
omipList.delete_experimentsByName_fromBuildersWithoutFlag(omipFileList, "-with-openmp")
omipList.print_list()


oceanList   = buildbot_experimentList("ocean")

myexps      = oceanList.add_experimentsByNameToBuildersWithoutFlag(\
[  "checksuite.ocean_internal/omip/exp.ocean_omip_ptest", \
   "checksuite.ocean_internal/technical/exp.test_ocean_omip_technical"], "without-mpi")
myexps      = oceanList.add_experimentsByNameToAllBuilders(["checksuite.ocean_internal/omip/exp.ocean_omip_short"])
myexps      = oceanList.add_experimentsByNameToMachine(["checksuite.ocean_internal/omip/exp.ocean_omip_testbed"], "mistral")

oceanList.print_list()

oceanList.write()

# some low level tests
sameOceanList  = buildbot_experimentList("ocean")
sameOceanList.read()
sameOceanList.print_list()

builderExps = sameOceanList.getBuilderExperimentNames("mistral_intel_hybrid")
print(builderExps)
quit()


myexp = sameOceanList.get_experiment("checksuite.ocean_internal/omip/exp.ocean_omip_short")
myexp.print_builders()
sameOceanList.delete_experiment(myexp)
# this will cause an error, since this experiment is deleted
myexp = sameOceanList.get_experiment("checksuite.ocean_internal/omip/exp.ocean_omip_short")


