#==============================================================================
# usage for interactve plot
#   SHOW=1 oceDiag.gp
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
set grid
set multiplot layout 2,2 title 'ICON OCEAN diagnostic'
plot 'oce_diagnostics.txt' using "step":"pot_energy" w l
set y2tics in
unset ytics
set grid y2tics
plot 'oce_diagnostics.txt' using "step":"kin_energy" w l axes x1y2
set ytics in
unset y2tics
set logscale y
plot 'oce_diagnostics.txt' using "step":"absolute_vertical_velocity" w l
set y2tics in
unset ytics
set grid y2tics
plot 'oce_diagnostics.txt' using "step":"total_salinity" w l axes x1y2, '' using "step":"total_temperature" w l axes x1y2
unset multiplot

# vim:ft=gnuplot
