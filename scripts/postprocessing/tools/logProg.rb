#!/usr/bin/env ruby
require 'pp'
require 'extcsv'
require 'jobqueue'
require 'gnuplot'
require 'date'

#===============================================================================
# usage:
#   OFILE='foo' OTYPE=pdf YMAX=100 Y2MAX=30 DEBUG=1 TITLE='title' logProg.rb <logFile list>
#
# All environment settings are optional:
#  OFILE: output file name without extension
#  OTYPE: image type, can be anythingm which is supported by gnuplot (e.g. png,svg,...)
#  YMAX: limits the first y-axes to given value (simulated years per day)
#  Y2MAX: limits the sec. y-axes to given value (read seconds per simulated day)
#  DEBUG: pring debug messages
#  TITLE: image title
#===============================================================================
#
# debug output
def dbg(msg); pp msg unless ENV['DEBUG'].nil? ; end
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
  dbg(simulationLengthInYears)

  # collect all requested timers
  dbg(logContent.grep(logpattern).map(&:split).map {|v|[v.grep(tagpattern),v[-1]].flatten})
  logContent.grep(logpattern).map(&:split).map {|v| [v.grep(tagpattern),v[-1]].flatten}.each {|k,v| 
    (retval[k] ||= [])<< v
    dbg(k)
  }
  # add numer of PEs
  retval['PEs'] = logContent.grep(/PE( |:)/).collect {|v| n = /PE(:| ) *(\d+)/.match(v)[2]}.uniq.size
  dbg(retval['PEs'])

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
def createPlot(data, oType='x11',oName='test')
  yearsPerDayData = ['total']
  # multiline plot with gnuplot
  Gnuplot.open do |gp|
    Gnuplot::Plot.new( gp ) do |plot|
      unless 'x11' == oType
        plot.terminal oType
        plot.output "#{oName}.#{oType}"
      end
      if ENV['TITLE'].nil?
        plot.title 'ICON ocean speed/runtime/scaling'
      else
        plot.title ENV['TITLE']
      end
      plot.grid

      plot.yrange  "[0:#{ENV['YMAX']}]"  unless ENV['YMAX'].nil?
      plot.y2range "[0:#{ENV['Y2MAX']}]" unless ENV['Y2MAX'].nil?

      plot.y2tics  'in'
      plot.key     'out horiz bot'
      plot.ylabel  'speed [years/day]'
      plot.y2label 'speed [real seconds/ simulated day]'
      data.each {|k,v|
        if yearsPerDayData.include?(k)
          plot.data << Gnuplot::DataSet.new( v ) do |ds|
            ds.with = "lines"
            ds.title = "#{k} (mean) [y/d]"
          end
          plot.data << Gnuplot::DataSet.new( v ) do |ds|
            ds.with = "errorbars"
            ds.title = "#{k} [y/d]"
          end
        else
          plot.data << Gnuplot::DataSet.new( v ) do |ds|
            ds.with = "errorbars axes x1y2"
            ds.title = "#{k} [s/d]"
          end
        end
      }
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

# Which timers should get analyzed?
tags    = %w[total gmres tracer_ab upd_phys upd_flx adv_horz dif_horz adv_vert dif_vert hflx_lim]
#tags    = %w[total global_sum exch_data]

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

data = {}
tags.each {|tag|
  data[tag] = a.datasets('PEs', tag, tag+'_min', tag+'_max').sort_by {|v| v[0].abs}.transpose
}

if ENV['OTYPE'].nil?
  createPlot(data)
else
  if ENV['OFILE'].nil?
    createPlot(data,ENV['OTYPE'])
  else
    createPlot(data,ENV['OTYPE'],ENV['OFILE'])
  end
end

