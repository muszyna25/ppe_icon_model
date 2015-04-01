# Author: ralf.mueller

$:.unshift File.join(File.dirname(__FILE__),".")
require 'list'
require 'minitest/autorun'
require 'pp'

OPERATORS = %w[mean min max div]
INTERVALS = %w[PT06H P01D P01M P01Y]
VARNAMES  = %w[u v temp zmc tend_t_up]

class MyListTest < Minitest::Test
  def setup
    @list = MyList.new
    @input = []
    OPERATORS.each {|operator| INTERVALS.each {|interval| VARNAMES.each {|varname|
          next if varname.size == 1 and interval == 'P01D'
          next if varname.size >  5 and interval == 'P01Y'
          next if varname.size == 4 and ['min','max'].include?(operator)
          @input << [operator,interval,varname]
    } } }
    #pp @input
  end

  def intervalHasEvent?(interval)
    return interval == (ENV.has_key?('INT') ? ENV['INT'] : 'P01D')
  end
  def loopA(output)
    output.each {|interval,joblist|
      if intervalHasEvent?(interval) then
        puts "INTERVAL:#{interval} is active"
        joblist.each {|operator,varlist|
          varlist.each {|var|
            puts "Work on VARNAME:#{var} for OPERATOR:#{operator}"
          }
        }
      end
    }
  end

  # multiple adds of the same value, should lead to a singe entry only
  def test_add_multiple
    @list.clear
    @list.add('a'); @list.add('a'); @list.add('a')
    assert_equal(1,@list.size)
    @list.add('b'); @list.add('b')
    assert_equal(2,@list.size)
  end

  # (A)
  # OUTPUT = {
  # "intervalA" => [
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vA,vC,...]},
  #    ],
  # "intervalB" => [
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vD,vC,...]},
  #    ],
  # }
  def test_output_A_hash
    output = {}
    @input.each {|line|
      operator,interval,varname = line
      ((output[interval] ||= {})[operator] ||= []) << varname
    }
    loopA(output)
  end


  # (B)
  # OUTPUT = {
  # "intervalA,operatorA" => [vA,vB,...],
  # "intervalA,operatorB" => [vA,vC,...],
  # "intervalB,operatorA" => [vA,vB,...],
  # "intervalB,operatorB" => [vD,vC,...],
  # }
  def test_output_B
  end
end
