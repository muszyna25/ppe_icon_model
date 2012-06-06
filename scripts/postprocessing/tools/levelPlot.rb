#!/usr/bin/env ruby
require 'extcsv'
require 'cdo'

# Usage
#
# ./levelPlot.rb <inputICONfile> <variableName> <operation:optional> <plotfile:optional>
#
#
# operations can be operators file fldmin (default), fldmax, fldavg, ... (everything that has only one input and one output stream)

ifile      = ARGV[0]
varname    = ARGV[1]
operation  = ARGV[2].nil? ? 'fldmin' : ARGV[2]
plotfile   = ARGV[3]
createPlot = (not ARGV[3].nil?)

#_plotfile  = "plot_#{varname}-#{Time.new.strftime("%Y%m%d-%H%M%S")}"

Cdo.debug = true

# Temporal file for text output
dataFile = MyTempfile.path

IO.popen("echo 'date|time|depth|#{varname}|' > #{dataFile}")
Cdo.outputkey('date,time,level,value', 
              :in => "-#{operation} -selname,#{varname} #{ifile} | sed -e 's/ 0 / 00:00:00 /' " +
                                                                "| sed -e 's/ 180000 / 18:00:00 /' " +
                                                                "| sed -e 's/ 120000 / 12:00:00 /' " +
                                                                "| sed -e 's/ 60000 / 06:00:00 /' | sed -e 's/ \\+/|/g' >>#{dataFile}")

icon = ExtCsv.new("file","psv",dataFile)

# Create datetime column for timeseries plot
icon.datetime = []
icon.date.each_with_index{|date,i| icon.datetime << [date,icon.time[i]].join(' ') }

pp icon.datetime

# Plot data with automatic splitting by depth
ExtCsvDiagram.plot_xy(icon,"datetime",varname.downcase,
                      "ICON: #{operation} on #{varname}", # Change title here
                      :label_position => 'below',:skipColumnCheck => true,
                      :type => 'lines',:groupBy => ["depth"], :onlyGroupTitle => true,
#                     :addSettings => ["logscale y"],     # Commend theses out for large scale values like Vert_Mixing_V
#                     :yrange => '[0.0001:10]',           # Otherwise you'll see nothing reasonable
                      :terminal => createPlot ? 'png' : 'x11',
                      :ylabel => "#{varname} [degC]",     # Correct the label if necessary
                      :input_time_format => "'%Y%m%d %H:%M:%S'",
                      :filename => plotfile,
                      :output_time_format => "'%d.%m.%y'",:size => "1600,600")
