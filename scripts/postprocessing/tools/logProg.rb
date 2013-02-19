#!/usr/bin/env ruby
require 'pp'
require 'extcsv'
require 'jobqueue'
require 'gnuplot'


#
# usage: 
#  procLog.rb logFile0 [ ... logFileN]

logFiles = ARGV

# preCheck files
logFiles.each {|file|
  unless File.exist?(file)
    warn "Cannot read file '#{file}'!"
    exit 1
  end
}
puts "Processing files #{logFiles.join(', ')} ..."

tags    = %w[total gmres tracer_ab upd_phys upd_flx]
q, lock = JobQueue.new, Mutex.new
results = {}

# setup processing method
procLog = lambda {|logFile,tags|
  logpattern = / (#{tags.join('|')}) /
  tagpattern = /(#{tags.join('|')})/
  retval     = {}

  # read in the whole file
  logContent = File.open(logFile).readlines.map(&:chomp)

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
    retval[k] = v_mean
    minMaxHash[k+'_min'] = v.min
    minMaxHash[k+'_max'] = v.max
  }
  retval.merge!(minMaxHash)

  #postProc, i.e. use arrays instead of numbers - this eases later plotting
  retval.each {|k,v| retval[k] = [v]}

  return retval
}

# process log files in parallel throuth the JobQueue
logFiles.each {|file|
  q.push {
    log = procLog[file,tags]
    lock.synchronize { results[file] = log }
  }
}
q.run

# collect all logs into a single object
a = ExtCsv.concat(*results.values.collect {|data| ExtCsv.new("hash","txt",data) })

##possible output for manual gnuplot visualization
#File.open("logPlot.dat","w") {|f| data4Plot.each {|v| f << v.join(' ') << "\n" } }

dataTotal = a.datasets('PEs', 'total', 'total_min', 'total_max').transpose
dataGmres = a.datasets('PEs', 'gmres', 'gmres_min', 'gmres_max').transpose
dataTrace = a.datasets('PEs', 'tracer_ab', 'tracer_ab_min', 'tracer_ab_max').transpose
# multiline plot with gnuplot
Gnuplot.open do |gp|
  Gnuplot::Plot.new( gp ) do |plot|
    plot.title 'ICON ocean runtime + scaling (181 days simulation time)'
    plot.grid
    plot.y2tics 'in'
    plot.data << Gnuplot::DataSet.new( dataTotal ) do |ds|
      ds.with = "errorbars"
      ds.title = 'total'
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
