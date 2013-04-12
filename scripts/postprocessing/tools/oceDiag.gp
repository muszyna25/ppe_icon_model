#==============================================================================
# usage for interactve plot
#   SHOW=1 oceDiag.gp
# user defined intut file
#   FILE=<ifile.txt> oceDiag.gp
#==============================================================================
#
# display the plot interactively, if SHOW is given on the command line
show=system("echo $SHOW")
if ( '' ne show ) set terminal x11 persist; else print(show); set terminal png size 1200,600; set output 'diag.png';

# get the input filename from the environment
file=system("echo $FILE")
if ( '' ne file ) print("use input file "); print(file); else file="oce_diagnostics.txt";

set grid
#set xdata time
#set timefmt "%Y-%m-%dT%H:%M%SZ"
#step datetime volume kin_energy pot_energy total_energy vorticity enstrophy potential_enstrophy absolute_vertical_velocity total_temperature total_salinity
set multiplot layout 3,2 title "ICON OCEAN diagnostic (".file.")"
plot file using 1:5 t "pot_energy" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 1:4 t "kin_energy" w l axes x1y2
set ytics in
unset y2tics
plot file using 1:10 t "absolute_vertical_velocity" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 1:12 t "total_salinity" w l
set logscale y
set ytics in
unset y2tics
plot file using 1:10 t "absolute_vertical_velocity" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 1:11 t "total_temperature" w l axes x1y2
unset multiplot

# vim:ft=gnuplot
