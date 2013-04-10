set terminal png size 1200,600
set output 'diag.png'
set grid
set multiplot layout 2,2 title 'ICON OCEAN diagnostic'
plot 'oce_diagnostics.txt' using "step":"pot_energy" w l
plot 'oce_diagnostics.txt' using "step":"kin_energy" w l
plot 'oce_diagnostics.txt' using "step":"absolute_vertical_velocity" w l
set logscale y
plot 'oce_diagnostics.txt' using "step":"total_salinity" w l, '' using "step":"total_temp" w l axes x1y2
unset multiplot
