$:.unshift File.join(File.dirname(__FILE__),"..","")
require 'test/unit'
require 'experiments'
require 'tempfile'
require 'pp'

class TestExperiments < Test::Unit::TestCase
  def test_nmlCreate
    nml = NmlLib.create('nml',:itopo => 0)
    assert_equal('nml',nml.keys.first)
    assert_equal({:itopo => 0},nml['nml'])
  end
  def test_nmlWrite
    nml      = lambda {|tag| NmlLib.create('nml_'+tag,:itopo => 0,:dz_lev => [50,50,50,100,200].map(&:to_f))}
    filename = Tempfile.new(self.class.to_s).path
    4.times {|i| NmlLib.write(nml.call(i.to_s),filename) }
    puts IO.popen("cat #{filename}").readlines
  end
end
