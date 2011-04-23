#! /bin/bash
#===========================================================================
# Description: create the grid image for the ICON page at the MPI-M website. 
# Author: Hui Wan (MPI-M, 2011-02-17)
#===========================================================================

 resolution="iconR2B04"
 grid_type="_spr0.90"
 extract_grid_data=0    # 1=YES,0=NO 

#===================================================================
# Generat grid files of GMT format from the standard ICON grid file
#===================================================================
 if [ $extract_grid_data -eq 1 ] ; then

    grid_data_path="/scratch/work/mh0287/hwan/exps/icon-dev/grids" 

    cdo outputbounds -selname,cell_area \
        ${grid_data_path}/${resolution}-grid${grid_type}.nc \
        > ${resolution}-cell_area.gmt
    cdo outputbounds -selname,dual_area \
        ${grid_data_path}/${resolution}-grid${grid_type}.nc \
        > ${resolution}-dual_area.gmt
 fi

#=============
# Plotting
#=============
 j=1
#for clat in -5 15 30 45 ; do    # try different center-latitutes
 for clat in 15 ; do 

   lonlat="-R-180/180/-90/90"    # plot the globe
   mapproj="-JG30/"$clat"/6i"    # use orthorgraphic projection
   output=landfill_${resolution}_$j

   # draw coast line, fill the land area; draw triangular and hexagonal grids

   pscoast $lonlat $mapproj -G186/208/196 -P -K >$output.ps
   psxy ${resolution}-cell_area.gmt -R -J -W0.2p,50/50/10 -P -M -O -K >>$output.ps
   psxy ${resolution}-dual_area.gmt -R -J -W0.1p,80/10/80 -P -M -O    >>$output.ps

  ## convert to pdf format
  #
  #convert $output.ps $output.pdf
  #rm $output.ps

   j=`expr $j + 1`
 done

#=============
