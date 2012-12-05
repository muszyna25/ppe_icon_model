#!/usr/bin/env ruby

require 'cdo'
require 'cdp'
require 'fileutils'
require 'jobqueue'
require 'socket'
require 'iconPlot'

#==============================================================================
# some helper methods =========================================================
def myPlotter
  plotFile = 'thingol' == Socket.gethostname \
           ? '/home/ram/src/git/icon/scripts/postprocessing/tools/icon_plot.ncl' \
           : '/pool/data/ICON/tools/icon_plot.ncl'
  plotter  = 'thingol' == Socket.gethostname \
           ? IconPlot.new(ENV['HOME']+'/local/bin/nclsh', plotFile, File.dirname(plotFile),'png','qiv',true,true) \
           : IconPlot.new("/home/zmaw/m300064/local/bin/nclsh", plotFile, File.dirname(plotFile), 'png','display',true,true)
  [plotter,plotFile]
end
#==============================================================================
def initFilename(experiment)
  "initial_#{experiment}.nc"
end
#==============================================================================
def secPlot(ofile,experiment,secPlots,q,lock,plotDir=".")
  plotDir << '/' unless '/' == plotDir[-1]

  title = (true) ? experiment : '"ICON Ocean, Mimetic-Miura, L40"'
  
  plotter, plotFile = myPlotter

  q.push {
    im = plotter.scalarPlot(ofile,plotDir+'T_'+     File.basename(ofile,'.nc'),'T',
                            :tStrg => title, :bStrg => '" "',
                            :hov => true,
                            :minVar => -3.0,:maxVar => 3.0,:withLines => false,:lStrg => 'T',
                            :numLevs => 20,:rStrg => 'Temperature', :colormap => "BlWhRe")
    lock.synchronize {(secPlots[experiment] ||= []) << im }
  }
  q.push {
    im =  plotter.scalarPlot(ofile,plotDir +'S_'+     File.basename(ofile,'.nc'),'S',
                             :tStrg => title, :bStrg => '" "',
                             :hov => true,
                             :minVar => -0.2,:maxVar => 0.2,:withLines => false,:lStrg => 'S',
                             :numLevs => 16,:rStrg => 'Salinity', :colormap => "BlWhRe")
    lock.synchronize {(secPlots[experiment] ||= []) << im }
  }
  q.push {
    im = plotter.scalarPlot(ofile,plotDir+'rhopot_'+File.basename(ofile,'.nc'),'rhopot',
                            :tStrg => title, :bStrg => '"  "',
                            :hov => true,
                            :minVar => -0.6,:maxVar => 0.6,:withLines => false,:lStrg => 'rhopot',
                            :numLevs => 24,:rStrg => 'Pot.Density', :colormap => "BlWhRe")
    lock.synchronize {(secPlots[experiment] ||= []) << im }
  }
end
#==============================================================================
def horizPlot(ofile,experiment,lastFile,plots,lock,plotDir=".")
  plotter, plotFile = myPlotter
  # compute the index of the last timestep
  lastTimestep = Cdo.ntime(:input => "-selname,ELEV #{lastFile}")[0].to_i - 1
  lastTimestepData = Cdo.seltimestep(lastTimestep, :input => "-selname,T #{lastFile}",:output => "lastTimeStep_"+File.basename(lastFile),:force => true)
  diffOfLast2Init = Cdo.sub(:input => [lastTimestepData,"-selname,T " +initFilename(experiment)].join(' '),:output => "diffOfLastTimestep_#{experiment}.nc")
    im = plotter.scalarPlot(diffOfLast2Init,'T_200m'+     File.basename(ofile,'.nc'),'T',
                            :tStrg => experiment, :bStrg => '" "',:maskName => "wet_c",:maskFile => lastFile,
                            :levIndex => 6,
                            :rStrg => 'Temperature')
    lock.synchronize {(plots[experiment] ||= []) << im }
    im = plotter.scalarPlot(diffOfLast2Init,'T_1000m'+     File.basename(ofile,'.nc'),'T',
                            :tStrg => experiment, :bStrg => '" "',:maskName => "wet_c",:maskFile => lastFile,
                            :levIndex => 12,
                            :rStrg => 'Temperature')
    lock.synchronize {(plots[experiment] ||= []) << im }
end
#==============================================================================
def cropSecPlots(secPlots,plotDir='.')
  plotDir << '/' unless '/' == plotDir[-1]
  q = JobQueue.new
  secPlots.each {|exp,files|
    q.push {
      cropfiles = []
      files.each {|sp|
        cropfile = plotDir+"crop_#{File.basename(sp)}"
        cmd      = "convert -resize 60%  -crop 650x560+50+50 #{sp} #{cropfile}"
        puts cmd
        system(cmd)
        cropfiles << cropfile
      }
      images = cropfiles[-3,3]
      system("convert +append #{(images.grep(/crop_T/) + images.grep(/crop_S/) + images.grep(/crop_rho/)).join(' ')} #{plotDir}#{exp}.png")
    }
  }
  q.run
  #system("convert +append #{(cropfiles.grep(/crop_T/) + cropfiles.grep(/crop_S/) + cropfiles.grep(/crop_rho/)).join(' ')} exp.png")
  #system("display #{cropfiles.join(' ')}") #if 'thingol' == Socket.gethostname
end
#==============================================================================
def cropMapPlots(plots,plotDir='.')
  plotDir << '/' unless '/' == plotDir[-1]
  q = JobQueue.new
  plots.each {|exp,files|
    q.push {
      cropfiles = []
      files.each {|sp|
        cropfile = plotDir+"crop_#{File.basename(sp)}"
        cmd      = "convert -crop 820x450+80+140 -resize 80% #{sp} #{cropfile}"
        puts cmd
        system(cmd)
        cropfiles << cropfile
      }
#      images = cropfiles[-3,3]
#      system("convert +append #{(images.grep(/crop_T/) + images.grep(/crop_S/) + images.grep(/crop_rho/)).join(' ')} #{plotDir}#{exp}.png")
    }
  }
  q.run
  #system("convert +append #{(cropfiles.grep(/crop_T/) + cropfiles.grep(/crop_S/) + cropfiles.grep(/crop_rho/)).join(' ')} exp.png")
  #system("display #{cropfiles.join(' ')}") #if 'thingol' == Socket.gethostname
end
#==============================================================================
#==============================================================================
# check input
if ARGV[0].nil?
  warn "no input files given"
  exit(-1)
end

files     = ( ARGV.size > 1 ) ? ARGV : Dir.glob(ARGV[0])
maskFile  = ENV["MASK"].nil? ? "mask.nc" : ENV["MASK"]
# check files
files.each {|file|
  warn "Cannot read file '#{file}'" unless File.exist?(file)
}
unless File.exist?(maskFile) 
  warn "Cannot open maskfile '#{maskFile}'"
  exit -1
end
#==============================================================================
q         = JobQueue.new([JobQueue.maxnumber_of_processors,16].min)
lock      = Mutex.new

Cdp.setCDO
Cdp.setDebug
Cdo.forceOutput = false
Cdo.debug       = true
diff2init       = true
def plot?
  not /(thingol|thunder)/.match(Socket.gethostname).nil?
end
#==============================================================================

# compute the experiments from the data directories and link the corresponding files
gridfile, experimentFiles, experimentAnalyzedData = Cdp.splitFilesIntoExperiments(files)

# process the files
#   start with selectiong the initial values from the first timestep
experimentFiles.each {|experiment, files|
  q.push {
    initFile    = initFilename(experiment)
    puts "Computing initial value file: #{initFile}"
    # create a separate File with the initial values
    if not File.exist?(initFile) or not Cdo.showname(:input => initFile).flatten.first.split(' ').include?("rhopot")
      initTS     = Cdo.selname('T,S',:input => "-seltimestep,1 #{files[0]}",:options => '-r -f nc')
      initRhopot = Cdo.rhopot(0,:input => initTS)
      merged = Cdo.merge(:input => [initTS,initRhopot].join(' '))
      FileUtils.cp(merged,initFile)
    end
  }
}
q.run
# compute meaked weight
maskedAreaWeights = Cdp.maskedAreaWeights("cell_area",gridfile,"wet_c",maskFile,"maskedAeraWeightsFrom#{File.basename(maskFile)}")

#   process the experiments results
experimentFiles.each {|experiment, files|
  files.each {|file|
    q.push {
      maskedYMeanFile = "masked_#{File.basename(file)}"
      fldmeanFile     = "fldmean_#{File.basename(file)}"
      rhopotFile      = "rhopot_#{File.basename(file)}"
      mergedFile      = "T-S-rhopot_#{File.basename(file)}"
      diffFile        = "T-S-rhopot_diff2init_#{File.basename(file)}"
      initFile        = initFilename(experiment)

      Cdo.div(:input => " -selname,T,S #{file} #{maskFile}",:output => maskedYMeanFile)
      # compute rhopot
      Cdo.rhopot(0,:input => maskedYMeanFile,:output => rhopotFile)

      Cdo.merge(:input => [maskedYMeanFile,rhopotFile].join(' '), :output => mergedFile)
      Cdo.sub(:input => [mergedFile,initFile].join(' '),:output => diffFile)
      Cdo.fldsum(:input => "-mul #{diffFile} #{maskedAreaWeights}", :output => fldmeanFile,:options => '-r -f nc')
      lock.synchronize {experimentAnalyzedData[experiment] << fldmeanFile }
    }
  }
}
q.run

# merge all yearmean data (T,S,rhopot) into one file per experiment
q.clear
secPlots, mapPlots = {}, {}
experimentAnalyzedData.each {|experiment,files|
  tag      = diff2init ? 'diff2init' : ''
  ofile    = [experiment,'T-S-rhopot',tag].join('_') + '.nc'

  yearmean = "yearmean"
  ymfile   = [yearmean,ofile].join("_")

  FileUtils.rm(ofile) unless plot? if File.exist?(ofile)
  FileUtils.rm(ymfile) unless plot? if File.exist?(ymfile)
  Cdo.cat(:input => files.sort.join(' '), :output => ofile, :force => !plot?)
  Cdo.settunits('years',:input => "-yearmean #{ofile}", :output => ymfile,:force => !plot?)

  ofile = ymfile
  q.push { secPlot(ofile,experiment,secPlots,q,lock) if plot? }
  q.push { horizPlot(ofile,experiment,experimentFiles[experiment][-1],mapPlots,lock) if plot? }
}
q.run
cropSecPlots(secPlots) if plot?
cropMapPlots(mapPlots) if plot?

