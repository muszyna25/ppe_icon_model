require 'tempfile'

module Cdo
  State = {}
  @@CDO = ENV['CDO'].nil? ? '/usr/bin/cdo' : ENV['CDO']

  def Cdo.Debug=(value)
    State[:debug] = value
  end
  def Cdo.Debug
    State[:debug]
  end
  def Cdo.getOperators
    cmd       = @@CDO + ' 2>&1'
    help      = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}
    if 5 >= help.size
      warn "Operators could not get listed by running the CDO binary (#{$opts[:bin]})"
      exit
    else
      help[help.index("Operators:")+1].split
    end
  end
  def Cdo.method_missing(sym, *args, &block)
    # args is expected to look like [opt1,...,optN,:in => iStream,:out => oStream] where
    # iStream could be another CDO call (timmax(selname(Temp,U,V,ifile.nc))
    puts "Operator #{sym.to_s} is called" if State[:debug]
    if getOperators.include?(sym.to_s)
      opts = args.empty? ? [] : ',' + args.reject {|a| a.class == Hash}.join(',')
      io = args.find {|a| a.class == Hash}
      run(" -#{sym.to_s}#{opts} #{io[:in]} ",io[:out])
    else
      "not found"
    end
  end
  def run(cmd,ofile=nil)
    if ofile == nil
      ofile = Tempfile.new('ifs2icon')
      ofile = ofile.path
    end 
    cmd =  "#{@@CDO} #{cmd} #{ofile}"
    call(cmd)
    return ofile
  end
  def call(cmd)
    if (State[:debug])
      puts '# DEBUG ====================================================================='
      puts cmd
      puts '# DEBUG ====================================================================='
      puts IO.popen(cmd).read
    else
      system(cmd + ' 1>/dev/null 2>&1 ')
    end
  end
end
if __FILE__ == $0
  include Cdo
  Cdo.Debug = false
  ifile = '/home/ram/data/examples/T.jan.nc'
  puts Cdo.seltimestep(1,:in => Cdo.selname('T',:in =>ifile),:out => 'ofile.nc')
#  puts Cdo.seltimestep(1,2,3,4,:in => 'in.nc',:out => 'ofile.nc')
end

