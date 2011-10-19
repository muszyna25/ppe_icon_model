require 'tempfile'
require 'thread'
require 'pp'
require 'extcsv'
require 'optparse'
require 'json'

# ==============================================================================
# CDO calling mechnism
module Cdo
  State = {}
  @@CDO = ENV['CDO'].nil? ? '/usr/bin/cdo' : ENV['CDO']

  # Only operators with documentation are accessible vie the build-in help.
  # Other have to be added manually
  @@undocumentedOperators = %w[geopotheight pressure_fl pressure_hl]
  @@addOperators          = %w[boundaryLevels thicknessOfLevels]

  def Cdo.Debug=(value)
    State[:debug] = value
  end
  def Cdo.Debug
    State[:debug]
  end

  # test if @@CDO can be used
  def Cdo.checkCdo
    unless (File.exists?(@@CDO) and File.executable?(@@CDO))
      warn "Testing application #@@CDO is not available!"
      exit 1
    else
      puts "Using CDO: #@@CDO"
      puts IO.popen(@@CDO + " -V").readlines
    end
  end

  def Cdo.setCdo(cdo)
    puts "Will use #{cdo} instead of #@@CDO" if Cdo.Debug
    @@CDO = cdo
  end


  def Cdo.getOperators
    cmd       = @@CDO + ' 2>&1'
    help      = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}
    if 5 >= help.size
      warn "Operators could not get listed by running the CDO binary (#{@@CDO})"
      pp help if Cdo.Debug
      exit
    else
      help[help.index("Operators:")+1].split
    end
  end
  def Cdo.call(cmd)
    if (State[:debug])
      puts '# DEBUG ====================================================================='
      puts cmd
      puts '# DEBUG ====================================================================='
      puts IO.popen(cmd).read
    else
      system(cmd + ' 1>/dev/null 2>&1 ')
    end
  end
  def Cdo.run(cmd,ofile=nil,options='')
    cmd = "#{@@CDO} -O #{options} #{cmd} "
    case ofile
    when $stdout
      cmd << " 2>/dev/null"
      return IO.popen(cmd).read.split
    when nil
      ofile = Tempfile.new("Ifs2Icon").path
    end
    cmd << "#{ofile}"
    call(cmd)
    return ofile
  end

  # Call an operator chain without checking opeartors
  def Cdo.chainCall(chain,*args)
    io = args.find {|a| a.class == Hash}
    args.delete_if {|a| a.class == Hash}

    chain   = chain.strip
    firstOp = chain
    firstOp = chain[0...[chain.index(','),chain.index(' ')].min] unless chain.index(',').nil?
    firstOp = firstOp[1..-1] if firstOp[0] == '-'
    if /(info|show|griddes)/.match(firstOp)
      Cdo.run(" #{chain} #{io[:in]} ",$stdout)
    else
      opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
      Cdo.run(" #{chain}#{opts} #{io[:in]} ",io[:out],io[:options])
    end
  end

  def Cdo.boundaryLevels(args)
    ilevels         = Cdo.showlevel(:in => args[:in]).map(&:to_f)
    bound_levels    = Array.new(ilevels.size+1)
    bound_levels[0] = 0
    (1..ilevels.size).each {|i| 
      bound_levels[i] =bound_levels[i-1] + 2*(ilevels[i-1]-bound_levels[i-1])
    }
    bound_levels
  end

  def Cdo.thicknessOfLevels(args)
    bound_levels = Cdo.boundaryLevels(args)
    delta_levels    = []
    bound_levels.each_with_index {|v,i| 
      next if i == 0
      delta_levels << v - bound_levels[i-1]
    }
    delta_levels
  end

  def Cdo.method_missing(sym, *args, &block)
    # args is expected to look like [opt1,...,optN,:in => iStream,:out => oStream] where
    # iStream could be another CDO call (timmax(selname(Temp,U,V,ifile.nc))
    puts "Operator #{sym.to_s} is called" if State[:debug]
    if getOperators.include?(sym.to_s) or @@undocumentedOperators.include?(sym.to_s)
      io = args.find {|a| a.class == Hash}
      args.delete_if {|a| a.class == Hash}
      if /(info|show|griddes)/.match(sym)
        run(" -#{sym.to_s} #{io[:in]} ",$stdout)
    else
        opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
        run(" -#{sym.to_s}#{opts} #{io[:in]} ",io[:out],io[:options])
    end
    else
      warn "Operator #{sym.to_s} not found"
  end
end
end
module MyTempfile
  require 'tempfile'
  @@_tempfiles = []
  @@persistent_tempfiles = false
  def MyTempfile.path
    unless @persistent_tempfiles
      t = Tempfile.new(self.class.to_s)
      @@_tempfiles << t
      t.path
    else
      t = "_"+rand(10000000).to_s
      @@_tempfiles << t
      t
    end
  end
end

