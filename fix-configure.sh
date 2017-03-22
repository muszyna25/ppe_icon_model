#! /bin/bash
#
# Luis Kornblueh, MPI-M, 2017-03-22

sed -i 's/|-LANG:=\* | -LIST:\* | -LNO:\* | -link)/|-LANG:=* | -LIST:* | -LNO:* | -link | -ltcmalloc | -ltcmalloc_minimal)/' configure
