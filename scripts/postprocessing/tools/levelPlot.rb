#!/usr/bin/env ruby
require 'extcsv'
require 'cdo'
require 'shellwords'

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

#Cdo.debug = true
Cdo.setCdo('/home/ram/src/cdo/trunk/cdo/build/bin/cdo')
#Cdo.checkCdo

# Temporal file for text output
dataFile = MyTempfile.path

# read the date
IO.popen("echo 'date|time|depth|#{varname}' > #{dataFile}")
Cdo.outputkey('date,time,level,value', 
              :in => "-#{operation} -selname,#{varname} #{ifile} >>#{dataFile}")
unit = Cdo.showunit(:in => "-selname,#{varname} #{ifile}").first

# postprocessing for correct time values
data = []
File.open(dataFile).each_with_index {|line,lineIndex|
  next if line.chomp.empty?
  _t = line.chomp.gsub(/ +/,'|').split('|')
  if 0 == lineIndex then
    data << _t
    next
  end
  if "0" == _t[1] then
    _t[1] = '00:00:00'
  else
    time = _t[1].reverse
    timeStr = ''
    while time.size > 2 do
      timeStr << time[0,2] << ':'
      time = time[2..-1]
    end
    timeStr << time.ljust(2,'0') unless time.size == 0
    _t[1] = timeStr.reverse
  end
  data << _t
}
icon = ExtCsv.new("array","plain",data.transpose)

# Create datetime column for timeseries plot
icon.datetime = []
icon.date.each_with_index{|date,i| icon.datetime << [date,icon.time[i]].join(' ') }
icon.datacolumns << "datetime"

# Plot data with automatic splitting by depth
unless icon.datacolumns.include?(varname)
  warn "Variable cannot be found!"
  exit -1
end
ExtCsvDiagram.plot_xy(icon,"datetime",varname,
                      "ICON: #{operation} on #{varname} (file:#{ifile})", # Change title here
                      :label_position => 'below',:skipColumnCheck => true,
                      :type => 'lines',:groupBy => ["depth"], :onlyGroupTitle => true,
#                     :addSettings => ["logscale y"],     # Commend theses out for large scale values like Vert_Mixing_V
#                     :yrange => '[0.0001:10]',           # Otherwise you'll see nothing reasonable
                      :terminal => createPlot ? 'png' : 'x11',
                      :ylabel => "#{varname} [#{Shellwords.escape(unit)}]",     # Correct the label if necessary
                      :input_time_format => "'%Y%m%d %H:%M:%S'",
                      :filename => plotfile,
                      :output_time_format => '"%d.%m.%y \n %H:%M"',:size => "1000,600")
