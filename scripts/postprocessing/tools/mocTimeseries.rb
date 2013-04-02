#!/usr/bin/env ruby
require 'jobqueue'
require 'cdo'

#==============================================================================
# USAGE:
# ./mocTimeseries.rb '../oce_mpiom_r11698_3600_1.0E-17_30/MOC*Z'
# 
#==============================================================================

mocPattern = ARGV[0].nil? ? 'MOC*Z' : ARGV[0]
mocVar     = ARGV[1].nil? ? 'var778' : ARGV[1]
title      = ENV['TITLE'].nil? ? mocPattern : ENV['TITLE']
plotFile   = ENV['PLOTFILE']
q          = JobQueue.new
lock       = Mutex.new
Cdo.debug  = true

# postprocessing/expanding
mocPattern = [File.expand_path(File.dirname(mocPattern)),File.basename(mocPattern)].join(File::SEPARATOR)
# calculate the moc timeseries in 1000 depth

mocFiles = Dir.glob(mocPattern).sort

levels = Cdo.nlevel(:input => "-seltimestep,1 #{mocFiles[0]}")[0].to_i

# process each moc file in parallel and join them at the end
mocOutputFiles = []
mocFiles.each_with_index {|file,i|
  q.push {
    levelSelection = 20 == levels ? "-sellevel,1000" : "-intlevel,1000 -sellevel,900/1100"
    ofile = Cdo.fldmean(:input => " -selname,#{mocVar} -mulc,1.e-9 -sellonlatbox,0,1,40,60 #{levelSelection} -yearmean -setgrid,r1x180 #{file}",
                :options => '-f nc -r')
                lock.synchronize { mocOutputFiles << [i,ofile]}
  }
}
q.run
pp mocOutputFiles
# merge together
mocComplete = Cdo.outputkey('date,value',:input => " -cat #{mocOutputFiles.sort.transpose[1].join(' ')}")
File.open("moc.dat","w") {|f| f << mocComplete.join("\n")}
IO.popen("LD_LIBRARY_PATH=/usr/lib gnuplot -p <<EOF
#{plotFile.nil? ? "" : "set terminal png size 800x400;set output '#{plotFile}'"}
set timefmt x '%Y-%m-%d'
set xdata time
set grid
set format x '%Y'
set title '#{title}' font 'arial,18pt'
plot 'moc.dat' using 1:2 w l t 'AMOC, 1000m, 50N (interpol. from 40N-60N)'
EOF").read
