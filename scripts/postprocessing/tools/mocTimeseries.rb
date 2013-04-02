#!/usr/bin/env ruby

require 'jobqueue'
require 'cdo'

mocPattern = ARGV[0].nil? ? 'MOC*Z' : ARGV[0]
mocVar     = ARGV[1].nil? ? 'var778' : ARGV[1]
q          = JobQueue.new
lock       = Mutex.new


# postprocessing/expanding
mocPattern = [File.expand_path(File.dirname(mocPattern)),File.basename(mocPattern)].join(File::SEPARATOR)
# calculate the moc timeseries in 1000 depth

mocFiles = Dir.glob(mocPattern)

levels = Cdo.nlevel(:input => "-seltimestep,1 #{mocFiles[0]}")[0].to_i

# process each moc file in parallel and join them at the end
mocOutputFiles = []
mocFiles.each {|file|
  q.push {
    levelSelection = 20 == levels ? "-sellevel,1000" : "-intlevel,1000 -sellevel,900/1100"
    ofile = Cdo.fldmean(:input => " -selname,#{mocVar} -mulc,1.e-9 -sellonlatbox,0,1,40,60 #{levelSelection} -yearmean -setgrid,r1x180 #{file}",
                :options => '-f nc -r')
                lock.synchronize { mocOutputFiles << ofile}
  }
}
q.run
pp mocOutputFiles
#LD_LIBRARY_PATH=/usr/lib gnuplot -p <<EOF
#set timefmt x "%Y-%m-%d"
#set xdata time
#set grid
#set format x "%Y"
#plot 'moc.dat' using 1:2 w l t 'AMOC, 1000m, 50N (interpol. from 40N-60N)'
#EOF
