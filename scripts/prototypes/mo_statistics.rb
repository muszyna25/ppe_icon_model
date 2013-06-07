# This is a prototype of a statistics module. Target language is fortran 95
require 'narray'

class DataStatistics

  # There are multiple cases of initialization
  # * use mode (Fixnum) to define the shape of the array, where the data is
  #   stored - values, which should be added must have the same shape. If mode
  #   is missing, then 
  # * the target argument has to be used. It defines an external NArray, where
  #   the data is accumulated. In addition to these two options
  # * the source paramters can be used to define an external NArray, from which
  # the data is read during the addtition to the stats object
  def initialize(target: [],source: [])
    @accumulation = target
    @source       = source
    @max, @min, @mean, @median = [],[],[],[]
    @numberOfAccumulations = 0
  end
  attr_accessor :accumulation, :max, :min ,:mean, :median, :numberOfAccumulations, :source

  # add values (singe, list or field) to the accumulation
  def add(values=nil)
    if @source.empty? and values.nil?
      warn "Nothing to add!"
      raise ArgumentError,
            "Both @source attribute and values parameters are not given!"
    end
    if (not @source.empty? and not values.nil?) then
      warn "To much to add!"
      raise ArgumentError,
            "Both @source attribute and values parameters are given!"
    end

    toBeAdded = values.nil? ? @source[0] : values

    @accumulation << toBeAdded
    nextValue = @accumulation.inject(0,&:+)
    @accumulation.clear
    @accumulation << nextValue

    @max          << toBeAdded.max  ; @min     << toBeAdded.min
    @mean         << toBeAdded.mean ; @median  << toBeAdded.median

    @numberOfAccumulations += 1
  end
end
class DataStatisticsCollection
  def initialze
    @statistics =[]
  end
  def add(statistic)
    @statistics << statistic
  end
  def update
    @statistics.each {|statistic|
      statistic.add([])
    }
  end
end

