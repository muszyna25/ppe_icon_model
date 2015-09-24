#!/usr/bin/env ruby

require './experiments'
include TestCases

exp = Experiment.new('vertival_mixing',
                     TestCase::TC33,
                     'R2B04',
                     [20,30,50,100,100,200,200,300,500,1000,1000,1000])

exp.run('squall')
