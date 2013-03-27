#!/usr/bin/env python
from pylab import *
from cdo import *
from string import rjust
import sys,os

#=============================================================================== 
# input file (expected is ICON ocean model result)
ifile     = sys.argv[1]
# variable name to plot (default T)
if (2 == sys.argv.__len__()):
  varname = 'T'
else:
  varname = sys.argv[2]
# location of the x-section
lon=[0,1]
lat=[-20,20]
# image scaling factor
scale = 5
# display debug stuff the script is started with DEBUG=<some value>
debug = 'DEBUG' in os.environ
#=============================================================================== 

cdo       = Cdo()
cdo.debug = debug

_ifile = cdo.sellonlatbox('%i,%i,%i,%i'%(lat[0],lat[1],lon[0],lon[1]),
                          input=' -remapnn,r360x180 -selname,%s %s'%(varname,ifile),
                          options='-P 8')

var = cdo.readCdf(_ifile).variables["T"]

shape = var.shape

rcParams['figure.figsize']=[shape[3]/scale,shape[1]/scale]

for t in range(0,shape[0]):
  ofile = 'xsect_'+rjust(str(t),4,'0')+'.png'
  print(ofile)
  v=var[t,0:shape[1],0,:]
  subplots_adjust(left=.05,right=.95,bottom=.1,top=.85)
  title('test')
  imshow(v,interpolation="nearest")
  #show()
  savefig(ofile,bbox_inches='tight',dpi=200)

