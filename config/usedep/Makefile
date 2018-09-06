CFLAGS=-std=c++11 -Wall -Wextra -Wunused -g -O2

INCFILES   := $(shell ls ../../src/*/*inc)
SRCFILES   := $(shell ls ../../src/*/*f90) $(shell ls ../../src/*/*/*f90)
JSBACHFILE := $(shell ls ../../src/lnd_phy_jsbach/*/*f90)
ALLFILES   := $(SRCFILES)
FILELIST   := $(filter-out ../../src/serialization/%.f90 ../../src/atm_dyn_iconam/mo_initicon.f90 $(JSBACHFILE),$(ALLFILES))

CXX        := g++
DEFINITIONS:= "__ICON__ __NO_JSBACH__  __NO_ICON_TESTBED__ __NO_ICON_ATMO__ __NO_HAMOCC__ __NO_ICON_PS_RAD__ ${EXTRADEF}"
SEARCHROOT := "../../src/drivers/icon.f90"

all: usedep_parser ocean.dep ocean.files
	cat ocean.dep

talktome:
	@echo "TALKTOME:: ${DEFINITIONS}"

usedep_parser:
	ragel usedep_parser.rl -o usedep_parser.cpp
	${CXX} ${CFLAGS} usedep_parser.cpp  -o usedep_parser -lm

clean:
	rm -f *.dot *~ *.o usedep_parser  usedep_parser.cpp 

check:
	echo $(FILELIST)
list:
	@time ./usedep_parser $(SEARCHROOT) $(DEFINITIONS) "$(FILELIST)"

ocean.files: usedep_parser
	@./usedep_parser $(SEARCHROOT) $(DEFINITIONS) "$(FILELIST)" > ocean.files

ocean.dep: usedep_parser
	@./usedep_parser $(SEARCHROOT) $(DEFINITIONS) "$(FILELIST)" -m > ocean.dep
