#!/usr/bin/env ruby
require 'tempfile'

# Example:
# ./globSum-over-time.rb iconOutput.nc T_global_sum.png 'T' 
#
# REQUIREMENTS: ruby, cdo, gnuplot
#
# Input data needs a reasonable time axis!

iFile       = ARGV[0]                                               # input file
oFile       = ARGV[1]                                               # output image file
varname     = ARGV[2].nil? ? 'T' : ARGV[2]                          # variable to process
title       = ARGV[3].nil? ? "global enery/inital energy" : ARGV[3] # graph title
globalTitle = ARGV[4].nil? ? '' : ARGV[4]                           # plot title

gnuplotScriptFile = varname+'_global.gpl'

tag='globalTend'

tfile0 = Tempfile.new(tag).path
tfile1 = Tempfile.new(tag).path
tfile2 = Tempfile.new(tag).path

cmd=<<END
cdo -fldsum  -vertsum -selname,#{varname} #{iFile} #{tfile0};
cdo infov -div #{tfile0} -timmax #{tfile0} > #{tfile1};
cat #{tfile1} | sed -e "s/ \\+/ /g" | cut -d ' ' -f 4,11 > #{tfile2};
END
gnuplotScript=<<END
set xdata time ; set timefmt "%Y-%m-%d"
set format x "%Y-%m"
set grid
set title '#{globalTitle}'
plot '#{tfile2}' using 1:2 w l title '#{title}'
set terminal png large size 800,400; set output '#{oFile}'
replot
END

puts cmd
puts '#==========================================='
puts gnuplotScript
puts '#==========================================='

system(cmd)
# write local gnuplot file for (evtl.) manual change
File.open(gnuplotScriptFile,"w") {|f| f << gnuplotScript}
puts IO.popen("gnuplot -persist #{gnuplotScriptFile}").read
