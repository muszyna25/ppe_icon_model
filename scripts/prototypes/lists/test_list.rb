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
          @input << [operator,interval,varname]
    } } } 
  end
  # multiple adds of the same value, should lead to a singe entry only
  def test_add_multiple
    @list.clear
    @list.add('a')
    @list.add('a')
    @list.add('a')
    assert_equal(1,@list.size)
    @list.add('b')
    @list.add('b')
    assert_equal(2,@list.size)
  end

  # OUTPUT = {
  # "intervalA" => [
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vA,vC,...]},
  #    ],
  # "intervalB" => 
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vD,vC,...]},
  #    ],
  #    }
  def test_output_L1

    @input.each {|line|
      pp line
    }
  end


  # OUTPUT = {
  # "intervalA,operatorA" => [vA,vB,...],
  # "intervalA,operatorB" => [vA,vC,...],
  # "intervalB,operatorA" => [vA,vB,...],
  # "intervalB,operatorB" => [vD,vC,...],
  # }
  def test_output_L2
  end
end
