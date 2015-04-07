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
    @input = []
    OPERATORS.each {|operator| INTERVALS.each {|interval| VARNAMES.each {|varname|
          next if varname.size == 1 and interval == 'P01D'
          next if varname.size >  5 and interval == 'P01Y'
          next if varname.size == 4 and ['min','max'].include?(operator)
          @input << [operator,interval,varname]
    } } }
  end

  def intervalHasEvent?(interval)
    return interval == (ENV.has_key?('INT') ? ENV['INT'] : 'P01D')
  end
  def loopHash(output)
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
  def check(listObj)
    listObj.clear
    listObj.add('a').add('a').add('a')
    assert_equal(1,listObj.size)
    listObj.add('b').add('b')
    assert_equal(2,listObj.size)
    pp listObj
  end
  def test_add_multiple
    check(MyHash.new)
    check(MyVector.new)
  end

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
  def test_output_hash
    output = {}
    @input.each {|line|
      operator,interval,varname = line
      ((output[interval] ||= {})[operator] ||= []) << varname
    }
    loopHash(output)
  end

  # OUTPUT = MyList [
  # ("intervalA", MyList[
  #   ("operatorA, MyList[vA,vB,...]),
  #   ("operatorB, MyList[vA,vc,...])
  #   ],
  # ),
  # ("intervalB", MyList[
  #   ("operatorA, MyList[vA,vB,...]),
  #   ("operatorB, MyList[vD,vc,...])
  #   ],
  # ),
  # ]
  def test_output_list
    #pp @input
   # output = My.new
    @input.each {|line|
      operator,interval,varname = line
    }
  end
end
