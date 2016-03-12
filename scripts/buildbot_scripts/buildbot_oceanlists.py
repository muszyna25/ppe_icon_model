#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# create the ocean experiment lists
#==============================================================================
from buildbot_builders import *
from model_paths import *
  

myPaths = model_paths()

omipFileList = myPaths.get_experimentsNames_inPaths(["checksuite.ocean_internal/omip/exp.*", "checksuite.ocean_internal/technical/exp.*"])
#print(omipFileList)
oceanList  = buildbot_experimentList("omip")
#myexps    = oceanList.add_experimentsByNameToAllBuilders(omipFileList)
#myexps    = oceanList.add_experimentsByNameToBuildersByName(omipFileList, ["mistral_gcc", "thunder_nag"])
#myexps    = oceanList.add_experimentsByNameToBuildersWithOptions(omipFileList, ["mistral", "daint_cpu"], None, None)
#oceanList.print_list()
#oceanList.delete_experimentsByNameFromBuildersWithOptions(omipFileList, ["mistral", "daint_cpu"], None, None)
#oceanList.print_list()
#myexps    = oceanList.add_experimentsByNameToBuildersWithOptions(omipFileList, None, ["openmp"], None)
myexps    = oceanList.add_experimentsByNameToBuildersWithOptions(omipFileList, ["mistral"], ["openmp"], None)
oceanList.print_list()
oceanList.delete_experimentsByName(["checksuite.ocean_internal/omip/exp.ocean_omip_long","checksuite.ocean_internal/omip/exp.ocean_omip_testSurface"])
oceanList.print_list()
quit()

oceanList.delete_experimentsByName_fromMachineName(omipFileList, "daint_cpu")
oceanList.print_list()
oceanList.delete_experimentsByName_fromMachineName(omipFileList, "thunder")
oceanList.print_list()
oceanList.delete_experimentsByName_fromBuildersWithoutFlag(omipFileList, "-with-openmp")
oceanList.print_list()


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


