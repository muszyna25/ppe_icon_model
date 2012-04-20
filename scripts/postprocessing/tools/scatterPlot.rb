#!/usr/bin/env ruby
require 'extcsv'
require 'cdo'
require 'jobqueue'
require 'pp'

# Usage
#
# ./scatterPlot.rb <inputICONfile> <variableNames:e.g. T,S> <timestep:optional> <plotfile:optional>
#
ifile      = ARGV[0]
varnames   = ARGV[1].split(',')
timestep   = ARGV[2].nil? ? '1' : ARGV[2]
plotfile   = ARGV[3]
createPlot = (not ARGV[3].nil?)
plotfile  = "plot_#{varnames.join}-#{Time.new.strftime("%Y%m%d-%H%M%S")}"


if varnames.size > 2
  warn "Please provide only to variable names"
  exit -1
end

#check, if region variable is available
regionName = 'rregio_c'
hasRegion = Cdo.showname(:in => ifile).join.split.include?(regionName)

# queue for parallel processing
jq = JobQueue.new

# output streams for read in values
outs = []
2.times do outs << [] end
outs << [] if hasRegion

# Get the data from the input file
variables = hasRegion ? [varnames,regionName].flatten : varnames
variables.each_with_index {|varname,i|
  jq.push {
    out  = outs[i]
    out << Cdo.outputkey('value', :in => "-selname,#{varname} -seltimestep,#{timestep} #{ifile}")
  }
}
jq.run

# Create the object
hash = hasRegion ? { varnames[0].to_sym => outs[0][0],
                     varnames[1].to_sym => outs[1][0],
                     regionName.to_sym  => outs[2][0]*10 } : {varnames[0].to_sym => outX[0],
                                                           varnames[1].to_sym => outY[0]}
icon = ExtCsv.new("hash","plain",hash)

#remove salinity values smaller than 1
[:s, :sal, :salinity].each {|col|
  if icon.datacolumns.include?(col.to_s)
    icon = icon.selectBy(col => "> 1")
    break
  end
  if icon.datacolumns.include?(col.to_s.upcase)
    icon = icon.selectBy(col.to_s.upcase.to_sym => "> 1")
    break
  end
}

# avoid gnuplot debugging output
$VERBOSE=false

# Plot data with automatic splitting by depth
ExtCsvDiagram.plot_xy(icon,varnames.first,varnames.last,
                      "ICON: Scatterplot on #{varnames.join(' and ')}", # Change title here
                      :groupBy => [:rregio_c],
                      :label_position => 'below',:skipColumnCheck => true,
                      :type => 'points', :onlyGroupTitle => true,
                      :terminal => createPlot ? 'png' : 'x11',
                      :ylabel => "#{varnames[1]}",     # Correct the label if necessary
                      :xlabel => "#{varnames[0]}",     # Correct the label if necessary
                      #                      :xrange => '[30:38]',
                      :filename => plotfile,
                      :size => "800,600")
# Plot an image fore each region separately
(icon.send(regionName).uniq.each {|r|
  jq.push {
    ExtCsvDiagram.plot_xy(icon.selectBy(:rregio_c => r),varnames.first,varnames.last,
                          "ICON: Scatterplot on #{varnames.join(' and ')}", # Change title here
                          :groupBy => [],
                          :label_position => 'below',:skipColumnCheck => true,
                          :type => 'points', :onlyGroupTitle => true,
                          :terminal => createPlot ? 'png' : 'x11',
                          :ylabel => "#{varnames[1]}",     # Correct the label if necessary
                          :xlabel => "#{varnames[0]}",     # Correct the label if necessary
                          #                      :xrange => '[30:38]',
                          :filename => plotfile+"_region-#{r}",
                          :size => "800,600")
  }
} and jq.run) if hasRegion
