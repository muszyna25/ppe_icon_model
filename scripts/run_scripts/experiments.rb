require 'pp'

# Methods for generating Namelist files
module NmlLib
  INDENT = '    '
  def NmlLib.create(name,values)
    {name => values}
  end
  def NmlLib.write(nml,filename)
    File.open(filename,"a") {|f|
      nml.each {|k,v|
        f << "&#{nml.keys.first}\n"
          v.each {|kk,vv|
            f << INDENT
            f << if vv.kind_of?(Array)
              "#{[kk,vv.join(',')].join(' = ')}\n"
            else 
              "#{[kk,vv].join(' = ')}\n"
            end
          }
        f << '/' << "\n"
      }
    }
  end
end

# Basic structure of what a testcase is
class TestCase
  attr_accessor :initialization, :forcing, :physics, :algorithms
  def initialize(name,initialization,forcing,algorithms)
    @name           = name
    @initialization = initialization
    @forcing        = forcing
    @algorithms     = algorithms
    @namelists      = []
  end
end
module TestCases
  TC33 = TestCase.new('33',nil,nil,nil)
end

# Implements the main functionality for running an experiment:
# * 
class Experiment
  @@defaulCompiler = 'gcc'
  @@defaultMpi     = 'openmpi'
  def initialize(name,testcase,resolution,verticalLevels,nProcs=8,nThreads=4)
    @name       = name
    @testcase   = testcase
    @resolution = resolution
    @nProcs     = nProcs
    @nThreads   = nThreads
  end
  def run(host,compiler=@@defaulCompiler,mpi=@@defaultMpi)
    puts "# create exp dir"
    puts "# setup and checkinitial data"
    puts "# setup and check forcing date"
    puts "# write name lists"
    pp @testcase
    puts "# load modules"
    puts "# run the model"
  end

  def setExperimentDir
  end
  def setupInputData
  end
  def setupForcingData
  end
  def checkFile(filename)
    checkHorizontalResolution
    checkVerticalLevels
  end
  def checkHorizontalResolution
  end
  def checkVerticalLevels
  end
  def createNamelists
  end
  private :setupForcingData, :setupInputData, :setExperimentDir, :checkFile, :checkHorizontalResolution, :checkVerticalLevels
end
