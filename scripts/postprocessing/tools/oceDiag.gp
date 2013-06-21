#==============================================================================
# usage for interactve plot
#   SHOW=1 gnuplot oceDiag.gp
# user defined intut file
#   FILE=<ifile.txt> gnuplot oceDiag.gp
#
# requirements: gnuplot 4.6 or later
#==============================================================================
#
# display the plot interactively, if SHOW is given on the command line
show=system("echo $SHOW")
if ( '' ne show ) {
  set terminal x11 persist
} else {
  print(show)
  set terminal png size 1200,600
  set output 'diag.png'
}
# get the input filename from the environment
file=system("echo $FILE")
if ( '' ne file ) {
  print("use input file ")
  print(file)
} else {
  file="oce_diagnostics.txt"
}
set grid
set xdata time
set timefmt '%Y-%m-%d %H:%M:%S'
set format x "%d.%m"
set multiplot layout 5,2 title "ICON OCEAN diagnostic (".file.")"
plot file using 2:"pot_energy" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 2:"kin_energy" w l axes x1y2
set ytics in
unset y2tics
plot file using 2:"absolute_vertical_velocity" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 2:"drake_passage" w l  axes x1y2
#set logscale y
set ytics in
unset y2tics
plot file using 2:"total_energy" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 2:"denmark_strait" w l axes x1y2
set ytics in
unset y2tics
plot file using 2:"total_energy" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 2:"gibraltar" w l axes x1y2
set ytics in
unset y2tics
plot file using 2:"ice_volume_nh" w l
set y2tics in
unset ytics
set grid y2tics
plot file using 2:"ice_extent_nh" w l axes x1y2
unset multiplot

# vim:ft=gnuplot
