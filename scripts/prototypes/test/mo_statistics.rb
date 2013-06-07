$:.unshift File.join(File.dirname(__FILE__),"..",".")
require "test/unit"
require "mo_statistics"
require "pp"

class TestCodeParser < Test::Unit::TestCase
  def test_create_single_stat
    # check plain scalar contructor
    (1..4).each {|i|
      ds = DataStatistics.new
      assert_equal([],ds.accumulation)
    }
  end
  def test_stat_with_target
    # check contructor with target
    target = []
    ds = DataStatistics.new(target: target)
    assert_equal(target, ds.accumulation)

    ds.add(NArray[[1.0,1.0],[1.0,1.0]])
    assert_equal(NArray[[1.0,1.0],[1.0,1.0]],ds.accumulation[0])

    ds.add(NArray[[1.0,1.0],[1.0,1.0]])
    assert_equal(NArray[[2.0,2.0],[2.0,2.0]],ds.accumulation[0])

    assert_equal(target,ds.accumulation) 
  end
  def test_stat_with_source
    # check contructor with source
    source = [NArray[[1.0,2.0],[3.0,4.0]]]
    ds = DataStatistics.new(source: source)

    ds.add
    assert_equal(source, ds.accumulation)
    ds.add
    assert_equal(2*source[0], ds.accumulation[0])
    ds.add
    ds.add
    assert_equal(4*source[0], ds.accumulation[0])

    source[0] = -NArray[[1.0,2.0],[3.0,4.0]]
    ds.add
    assert_equal(-3*source[0], ds.accumulation[0])
  end
  def test_stat_with_source_and_target
    source = [NArray.float(2,2)+1]
    target = [NArray.float(2,2)+10]
    ds = DataStatistics.new(source: source, target: target)

    ds.add
    assert_equal([NArray.float(2,2)+11],ds.accumulation)

    source[0] = NArray.float(2,2)+2
    ds.add
    assert_equal([NArray.float(2,2)+13],ds.accumulation)
  end
  def test_init
    test_create_single_stat
    test_stat_with_target
    test_stat_with_source
    test_stat_with_source_and_target
  end
  def test_addTo_single_stat
    addMe = NArray[1.0,2.0,3.0]
    ds = DataStatistics.new
    ds.add(addMe)
    assert_equal(1,ds.numberOfAccumulations)
    ds.add(addMe)
    assert_equal(2,ds.numberOfAccumulations)
    10.times { ds.add(addMe) }
    assert_equal(12,ds.numberOfAccumulations)
    assert_equal(12.0,ds.accumulation.first[0])
    assert_equal(24.0,ds.accumulation.first[1])
    assert_equal(36.0,ds.accumulation.first[2])

    assert_equal(1.0,ds.min.uniq[0])
    assert_equal(2.0,ds.mean.uniq[0])
    assert_equal(3.0,ds.max.uniq[0])
    assert_equal(2.0,ds.median.uniq[0])
  end
  def test_add_usage
    ds = DataStatistics.new
    assert_raise ArgumentError do; ds.add; end
    ds = DataStatistics.new(source: [1])
    assert_raise ArgumentError do; ds.add(2); end
  end

  def create_stat_collections
  end
  def update_stat_collection
  end
end
