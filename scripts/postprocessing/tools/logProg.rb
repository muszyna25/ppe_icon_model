#!/usr/bin/env ruby
require 'pp'
require 'extcsv'
require 'jobqueue'
require 'gnuplot'
require 'date'

#===============================================================================
#
# proccesing method which should be called for each file
def procLog(logFile,tags)
  logpattern = / (#{tags.join('|')}) /
  tagpattern = /(#{tags.join('|')})/
  retval     = {}

  # read in the whole file
  logContent = File.open(logFile).readlines.map(&:chomp)

  # get the simulation length
  simulationStart, simulationEnd = logContent.grep(/(start|end)_date/).map {|v| DateTime.parse(v.split(' = ')[-1])}
  simulationLengthInDays         = (simulationEnd - simulationStart).to_i
  simulationLengthInYears        = simulationLengthInDays/365.0
  oneDayInSeconds                = 86400.0

  # collect all requested timers
  logContent.grep(logpattern).map(&:split).map {|v| [v[0,3].grep(tagpattern),v[-1]].flatten}.each {|k,v| 
    (retval[k] ||= [])<< v
  }
  # add numer of PEs
  retval['PEs'] = logContent.grep(/PE( |:)/).collect {|v| n = /PE(:| ) *(\d+)/.match(v)[2]}.size
  # use min,max,mean instead of all available values
  minMaxHash = {}
  retval.each {|k,v|
    next if 'PEs' == k
    v = v.map(&:to_f)
    v_mean = v.reduce(0.0,:+)/v.size
    if 'total' == k
      retval[k] = simulationLengthInYears/(v_mean/oneDayInSeconds)
      minMaxHash[k+'_min'] = simulationLengthInYears/(v.min / oneDayInSeconds)
      minMaxHash[k+'_max'] = simulationLengthInYears/(v.max / oneDayInSeconds)
    else
      retval[k] = v_mean / simulationLengthInDays
      minMaxHash[k+'_min'] = v.min / simulationLengthInDays
      minMaxHash[k+'_max'] = v.max / simulationLengthInDays
    end
  }
  retval.merge!(minMaxHash)

  #postProc, i.e. use arrays instead of numbers - this eases later plotting
  retval.each {|k,v| retval[k] = [v]}

  return retval
end

# create a plot for some dedicated parameters: total runtime, gmres and tracer advection
def createPlot(dataTotal,dataGmres,dataTrace,oType='x11',oName='test')
  # multiline plot with gnuplot
  Gnuplot.open do |gp|
    Gnuplot::Plot.new( gp ) do |plot|
      unless 'x11' == oType
        plot.terminal oType
        plot.output "#{oName}.#{oType}"
      end
      plot.title 'ICON ocean speed/runtime/scaling R2B04 on blizzard'
      plot.grid
      plot.y2tics  'in'
      plot.key     'out horiz bot'
      plot.ylabel  'speed [years/day]'
      plot.y2label 'speed [seconds/day]'
      plot.data << Gnuplot::DataSet.new( dataTotal ) do |ds|
        ds.with = "lines"
        ds.title = 'total'
      end
      plot.data << Gnuplot::DataSet.new( dataTotal ) do |ds|
        ds.with = "errorbars"
        ds.title = 'total error'
      end
      plot.data << Gnuplot::DataSet.new( dataGmres ) do |ds|
        ds.with = "errorbars axes x1y2"
        ds.title = 'gmres'
      end
      plot.data << Gnuplot::DataSet.new( dataTrace ) do |ds|
        ds.with = "errorbars axes x1y2"
        ds.title = 'tracer_ab'
      end
    end
  end
end

#===============================================================================
#====== MAIN PROGRAM ===========================================================
# usage: 
#  procLog.rb logFile0 [ ... logFileN]
#

logFiles = ARGV

# preCheck file existenc
logFiles.each {|file|
  unless File.exist?(file)
    warn "Cannot read file '#{file}'!"
    exit 1
  end
}
puts "Processing files #{logFiles.join(', ')} ..."

# Which timer should get analyzed?
tags    = %w[total gmres tracer_ab upd_phys upd_flx]

# setup for parallel processing
q, lock = JobQueue.new, Mutex.new
results = {}

# process log files in parallel throuth JobQueue
logFiles.each {|file|
  q.push {
    log = procLog(file,tags)
    lock.synchronize { results[file] = log }
  }
}
q.run

# collect all logs into a single object
a = ExtCsv.concat(*results.values.collect {|data| ExtCsv.new("hash","txt",data) })

##possible output for manual gnuplot visualization
#File.open("logPlot.dat","w") {|f| data4Plot.each {|v| f << v.join(' ') << "\n" } }

dataTotal = a.datasets('PEs', 'total', 'total_min', 'total_max').sort_by {|v| v[0].abs}.transpose
dataGmres = a.datasets('PEs', 'gmres', 'gmres_min', 'gmres_max').transpose
dataTrace = a.datasets('PEs', 'tracer_ab', 'tracer_ab_min', 'tracer_ab_max').transpose

createPlot(dataTotal,dataGmres,dataTrace)
createPlot(dataTotal,dataGmres,dataTrace,'png')

