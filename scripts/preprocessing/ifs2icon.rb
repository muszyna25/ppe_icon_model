#!/usr/bin/env ruby
# Due to the fact, that it will make heavy use of CDO, there is an CDO issue
# for that topic: https://code.zmaw.de/issues/612
#
# Use 'rdoc ifs2icon.rb' to get documentation

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
    puts "Will use #{cdo} instead of #@@CDO"
    @@CDO = cdo
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
  def Cdo.method_missing(sym, *args, &block)
    # args is expected to look like [opt1,...,optN,:in => iStream,:out => oStream] where
    # iStream could be another CDO call (timmax(selname(Temp,U,V,ifile.nc))
    puts "Operator #{sym.to_s} is called" if State[:debug]
    if getOperators.include?(sym.to_s)
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

  # Call an operator chain without checking opeartors
  def Cdo.chainCall(chain,*args)
    io = args.find {|a| a.class == Hash}
    args.delete_if {|a| a.class == Hash}
    if /(info|show|griddes)/.match(sym)
      run(" -#{chain} #{io[:in]} ",$stdout)
    else
      opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
      run(" -#{chain}#{opts} #{io[:in]} ",io[:out],io[:options])
    end
  end
end
# ==============================================================================
# Option handling
class PreProcOptions
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    options                    = OpenStruct.new
    options.interpolation_type = :simple
    options.verbose            = false
    options.openmp             = 2

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: ifs2icon.rb [options]"

      opts.separator ""
      opts.separator "Specific options:"

      # Mandatory arguments
      opts.on("-g", "--grid-file FILE",
              "Use FILE read grid definition/resolution") do |file|
        options.gridfile = file
              end
      opts.on("-i", "--input-file FILE",
              "Use FILE read input variables") do |file|
        options.inputfile = file
              end
      opts.on("-o", "--output-file FILE",
              "Use FILE write output variables") do |file|
        options.outputfile = file
              end

      # optional configuration file using csv like format
      opts.on("-c", "--config-file FILE", "Use FILE as configuration file", "Format should be a separated columns text file (e.g. csv)") do |file|
        options.configurationfile = file
      end
      # optional configuration file using JSON like format
      jsonStringExample=<<-EOF
                  [
                    ["Description"       , "inputname", "outputname", "code", "grid", "typeOfLayer", "nlevel", "GP", "I", "Notes"],
                    ["temperature "      , "T"        , "T"         , "130" , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["zonal wind comp. u", "U"        , "U"         , "131" , "cell", "hybridLayer", "nlev"  , "2" , "Q", ""]     ,
                    ["normal velocity"   , "VN "      , ""          , ""    , "edge", "hybridLayer", "nlev"  , "2" , "" , ""]     ,
                    ["vertical velocity" , "OMEGA"    , "W"         , "135" , "cell", "hybridLayer", "nlev+1", "1" , "" , ""]     ,
                    ["Pressure"          , "P"        , "PRES"      , ""    , "cell", "hybridLayer", "nlev"  , "1" , "" , ""]     ,
                    ["air density"       , "RHO"      , ""          , ""    , "cell", "hybridLayer", "nlev"  , "1" , "" , ""]     ,
                    ["exner pressure "   , "EXNER"    , ""          , ""    , "cell", "hybridLayer", "nlev"  , "1" , "" , ""]     ,
                    ["specific humidity" , "QV"       , "Q"         , "133" , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["cloud ice content" , "QI"       , "CIWC"      , "247" , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["rain water content", "QR"       , "CRWC"      , "75"  , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["snow water content", "QS"       , "CSWC"      , "76"  , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["ozone mixing ratio", "O3"       , "O3"        , "203" , "cell", "hybridLayer", "nlev"  , "1" , "Q", ""]     ,
                    ["surface pressure"  , "PS"       , "LNSP"      , "152" , "cell", "surface"    , "1"     , "1" , "Q", ""]     ,
                    ["snow temperature"  , "T_SNOW"   , "TSN"       , "238" , "cell", "surface"    , "1"     , "1" , "M", ""]     ,
                    ["density of snow"   , "RHO_SNOW" , "RSN"       , "33"  , "cell", "surface"    , "1"     , "1" , "M", ""]     ,
                    ["surface roughness" , "Z0"       , "SR"        , "173" , "cell", "surface"    , "1"     , "1" , "M", ""]
                  ]
                EOF
      opts.on("-j", "--json-config-file FILE", "Use FILE as JSON configuration file","Format example:\n"+jsonStringExample) do |data|
        options.configurationdata = data
      end

      # Optional argument with keyword completion.
      opts.on("-t", "--interpolation-type [TYPE]", [:simple, :bilinear, :bicubic],
              "Select interpolation type (simple, bilinear, bicubic)") do |t|
        options.interpolation_type = t
              end
      #
      # List of output variables
      opts.on("-O", "--output-variables x,y,z", Array, "'list' of output variables") do |list|
        options.outputVars = list
      end

      # Set OpenMP multithreadding for CDO
      opts.on("-P", "--openmp p", Numeric, "Set number of OpenMP threads to <p>") do |p|
        options.openmp = p
      end

      # Boolean switches
      opts.on("-s", "--strict", "Process only variables which have valid entry in the configuration") do |v|
        options.strict = v
      end
      opts.on("-C", "--check", "Check input file and configuration on which output variables are available.") do |v|
        options.check = v
      end
      opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
        options.verbose = v
      end
      opts.on("-D", "--debug", "Run in debug mode") do |v|
        options.debug = v
      end

      # create zlevel options for each given output variables
      options.zlevels = {}
      opts.on("-z","--zlevels varName0,lev0,lev1,levN,varName1,lev0,...",
              Array,
              "'list' of target zlevel values for the output variables") {|list|
                list.slice_before {|v| v.to_i.to_s != v}.each {|sublist|
                  options.zlevels[sublist[0]] = sublist[1..-1]
                }
      }

      opts.separator ""
      opts.separator "Common options:"

      # No argument, shows at tail.  This will print an options summary.
      # Try it and see!
      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end

      # Another typical switch to print the version.
      opts.on_tail("--version", "Show version") do
        puts 'PreProc ' + [0,0,1].join('.')
        exit
      end
    end

    opts.parse!(args)

    begin
      mandatory = [:outputfile, :inputfile, :gridfile]
      missing = mandatory.select{ |param| options.send(param).nil? }
      if not missing.empty?
        warn "Missing options: #{missing.join(', ')}"
        puts opts
        exit
      end
    rescue OptionParser::InvalidOption, OptionParser::MissingArgument
      puts $!.to_s
      puts opts
      exit
    end
    options
  end
end
# ==============================================================================
# variable substitution between icon and Eecmwf
module Ecmwf2Icon
  _defaultConfig =[
    ['Description'                        , 'outputname', 'inputname', 'code', 'grid', 'typeOfLayer', 'nlevel', 'GP', 'I', 'Notes']                                        ,
    ['temperature '                       , 'T'         , 'T'        , '130' , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['zonal wind comp. u'                 , 'U'         , 'U'        , '131' , 'cell', 'hybridLayer', 'nlev'  , '2' , 'Q', '']                                             ,
    ['meridional wind comp. v '           , 'V'         , 'V'        , '132' , 'cell', 'hybridLayer', 'nlev'  , '2' , 'Q', '']                                             ,
    ['normal velocity'                    , 'VN '       , ''         , ''    , 'edge', 'hybridLayer', 'nlev'  , '2' , '' , 'comp. using U and V']                          ,
    ['vertical velocity'                  , 'OMEGA'     , 'W'        , '135' , 'cell', 'hybridLayer', 'nlev+1', '1' , '' , 'hydr. Approximation Omega -> w (Pa/s -> m/s)r'],
    ['Pressure'                           , 'P'         , 'PRES'     , ''    , 'cell', 'hybridLayer', 'nlev'  , '1' , '' , '? use pressure or density ? if pressure        , for vertical interpolation take pressure deviation from reference pp=pres-ref'],
    ['air density'                        , 'RHO'       , ''         , ''    , 'cell', 'hybridLayer', 'nlev'  , '1' , '' , 'use gas eq. RHO=f(T                            , p                                                                              , (qv))'],
    ['virtual potential temperature'      , 'THETA_V'   , ''         , ''    , 'cell', 'hybridLayer', 'nlev'  , '1' , '' , 'to be computed inside ICON']                   ,
    ['exner pressure '                    , 'EXNER'     , ''         , ''    , 'cell', 'hybridLayer', 'nlev'  , '1' , '' , '']                                             ,
    ['specific humidity'                  , 'QV'        , 'Q'        , '133' , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['cloud liquid water content'         , 'QC'        , 'CLWC'     , '246' , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['cloud ice content'                  , 'QI'        , 'CIWC'     , '247' , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['rain water content'                 , 'QR'        , 'CRWC'     , '75'  , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['snow water content'                 , 'QS'        , 'CSWC'     , '76'  , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['ozone mixing ratio'                 , 'O3'        , 'O3'       , '203' , 'cell', 'hybridLayer', 'nlev'  , '1' , 'Q', '']                                             ,
    ['surface pressure'                   , 'PS'        , 'LNSP'     , '152' , 'cell', 'surface'    , '1'     , '1' , 'Q', '']                                             ,
    ['specific humidity at surface'       , 'QV_S'      , ''         , ''    , 'cell', 'surface'    , '1'     , '1' , 'M', 'take Q at lowest level ?']                     ,
    ['snow temperature'                   , 'T_SNOW'    , 'TSN'      , '238' , 'cell', 'surface'    , '1'     , '1' , 'M', '']                                             ,
    ['water content of snow'              , 'W_SNOW'    , 'SD'       , '141' , 'cell', 'surface'    , '1'     , '1' , 'M', '']                                             ,
    ['density of snow'                    , 'RHO_SNOW'  , 'RSN'      , '33'  , 'cell', 'surface'    , '1'     , '1' , 'M', '']                                             ,
    ['water cont. of interception storage', 'W_I'       , 'SRC'      , '198' , 'cell', 'surface'    , '1'     , '1' , 'M', '']                                             ,
    ['surface roughness'                  , 'Z0'        , 'SR'       , '173' , 'cell', 'surface'    , '1'     , '1' , 'M', '']
  ]
  DefaultConfig = ExtCsv.new('array','plain',_defaultConfig.transpose)

  DefaultColumnNames =
    {
      :outName                 => :outputname,
      :inName                  => :inputname,
      :code                    => :code,
      :grid                    => :grid,
      :levType                 => :typeOfLayer,
      :nlevels                 => :nlevel,
      :info                    => :Description,
      :horizontalInterpolation => :I
    }

  DefaultInputVars = {
    135 => nil, # W vertival velocity
    129 => nil, # Z geopotential
    130 => nil, # T Temperature
    133 => nil, # Q spec. humidity
    138 => nil, # vor vorticity
    152 => nil, # lnsp Logarithm of surface pressure
    155 => nil, # div divergence
    203 => nil  # o3 ozone mass mixing ratio
  }

  DefaultOutputVars = %w[U V T]

  State = {}

  def defaultConfig
    DefaultConfig
  end

  def readConfig(configFile,itype='txt')
    case itype
    when 'txt'
      config = ExtCsv.new("file","txt",configFile)
      unless Ecmwf2Icon.check(@config)
        warn "Cannot read in file #{configFile}! Please check its structure"
        exit
      end
    when 'json'
      config = ExtCsv.new('array','plain',JSON.parse(File.open(configFile,"r").read).transpose)
    else
      warn 'Unusable input type for configuration data'
      exit
    end
    config
  end

  def Ecmwf2Icon.displayConfig(config,verbose=false, debug=false)
    Dbg.msg(config.columns(DefaultColumnNames[:code],
                           DefaultColumnNames[:outName],
                           DefaultColumnNames[:inName],
                           DefaultColumnNames[:info]).to_string("tsv",false),
            verbose,
            debug)
  end

  def Ecmwf2Icon.hasCode?(code)
    return true if @config.send(DefaultColumnNames[:code]).include?(code)
    return false
  end

  def Ecmwf2Icon.grid(varid,conf)
    case
    when varid.kind_of?(Numeric)
      return conf.selectBy(DefaultColumnNames[:code] => "'#{varid}'").grid[0]
    when varid.kind_of?(String)
      return conf.selectBy(DefaultColumnNames[:outName] => varid).grid[0]
    end
    warn "Could not find grid information for var '#{varid}'!"
    warn "Check config with 'checkConfig'"
    exit
  end

  def Ecmwf2Icon.hasGridinfo?(varid,conf)
    case
    when varid.kind_of?(Numeric)
      return true if %w[cell edge vertex].include?(conf.selectBy(DefaultColumnNames[:code] => "'#{varid}'").grid[0])
    when varid.kind_of?(String)
      return true if %w[cell edge vertex].include?(conf.selectBy(DefaultColumnNames[:outName] => varid).grid[0])
    end
    return false
  end

  def Ecmwf2Icon.checkConfig(config)
    Ecmwf2Icon.checkGriddes(config)
  end

  def Ecmwf2Icon.checkGriddes(configObject)
    return false unless configObject.send(DefaultColumnNames[:grid]).uniq.sort == %w[cell edge vertex]
    return true
  end
# def checkNames(configObject)
#   # every codenummer must have valid names
#   codes = configObject.id.delete_if {|code| code.to_i == 0 and code != '0'}
# end
end
# ==============================================================================
# Debuging facilities
module Dbg
  def Dbg.msg(str,verboseMode=false,debugMode=false)
    if debugMode
      puts str
      return
    elsif verboseMode
      puts str
    end
  end
end

# ==============================================================================
# Sized Queue for limiting the number of parallel jobs
class JobQueue
  attr_reader :size, :queue, :threads

  def initialize(size)
    @size  = size
    @queue = Queue.new
  end

  def push(*items)
    items.flatten.each {|it| @queue << it}
  end

  def run
    @threads = (1..@size).map {|i|
      Thread.new(@queue) {|q|
        until ( q == ( task = q.deq ) )
          system(task)
        end
      }
    }
    @threads.size.times { @queue.enq @queue}
    @threads.each {|t| t.join}
  end
end
# example script
# if $0 == __FILE__
#   size = 8
#   cmds = (0...1000).map {|i|
#     "cdo remapnn,r10x10 /home/ram/data/examples/T.jan.nc ~/#{i}.nc"
#   }
#   q = JobQueue.new(size)
#   q.push(cmds)
#   q.run
# end

class PreProc; end
# ==============================================================================
# A preprocessor should basically constist fof the following parts
# * knowledge about required in and output variables
# * rules for computing the output from the input variables
# * abtraction from the IO level
# The folloing steps have to be done
# * read input variables
# * compute the output variables
# * map output variables onto the required grid
# * (optional) apply land-see-mask
# * (optional) to concistency checks

class Ifs2Icon

  include Ecmwf2Icon

  # These variables are dummies to make the 3 different grids accessable to
  # CDO. They are created by ICON's grid generator
  GRID_VARIABLES = {
    :cell   => :ifs2icon_cell_grid,
    :edge   => :ifs2icon_edge_grid,
    :vertex => :ifs2icon_vertex_grid
  }

  # Intermediate filenames are created multithreaded. This  prevents them  from
  # beeing garbage collected and enables a threadsafe use of ruby's tempfile
  # library
  @@_tempfiles = []

  attr_accessor :inVars, :outVars, :rules, :config, :ifile

  def initialize(options)

    @lock = Mutex.new

    # required computation for creating the output variables
    # TODO this is still a dummy for later implentation steps
    @rules    = {}

    @options = options
    @options.freeze

    @tempfiles = []

    @inVars, @outVars, @gridVars = {}, {}, {}

    # Create internal configuration state
    if @options.configurationfile
      Dbg.msg('Use configuration data from TXT file',@options.verbose)
      @config = readConfig(@options.configurationfile)
    elsif @options.configurationdata
      Dbg.msg('Use configuration data from JSON file',@options.verbose)
      @config = readConfig(@options.configurationdata,'json')
    else
      Dbg.msg('Use default configuration',@options.verbose)
      @config = Ecmwf2Icon::DefaultConfig
    end

    Ecmwf2Icon.checkConfig(@config)

    Ecmwf2Icon.displayConfig(@config,@options.verbose, @options.debug)

    prepareGridfiles

    @ifile = @options.inputfile

    # Read input variables form File
    @inVars = getInputVariables

    # Check which input variabels have a valid configuration
    checkConfigOnInputVariables
  end

  def checkForPossibleOutput
    ovars = checkConfigOnInputVariables
    puts "Possible output variables are #{ovars.join(',')}" 
  end

  def run
    readVars

    # Check output variables form the cmdline option if they have valid grid
    # info in the configuration
    @outVars = getAndCheckOutputVariables if @outVars.empty?

    computeOutputVars

    applyLSM

    consistencyCheck

    writeOutput(@options.outputfile)
  end

  # Create files for later usage as grid description for the remapping
  def prepareGridfiles
    GRID_VARIABLES.each {|k,v|
      Dbg.msg("Selecting #{k} grid from #{@options.gridfile}",false, @options.debug)
      varfile = tfile
      Cdo.selname(v,:in => @options.gridfile,:out => varfile)
      @gridVars[v] = varfile
      Dbg.msg("#{k} grid is now in #{varfile}",false, @options.debug)
      Dbg.msg(IO.popen("ls -crtlh #{varfile}").read,false,@options.debug)
    }
  end

  def getInputVariables
    inVars = {}
    # Use names for netcdf files and codes for grib input files
    operator  = /grb/.match(File.extname(@ifile)) ? :showcode : :showname
    inputVars = Cdo.send(operator,:in => @ifile)
    # Use numbers as identifiers for grib input
    inputVars.map!(& :to_i) if operator == :showcode

    Dbg.msg("Found these variables in the input file = #{inputVars.join(",")} in input file",@options.verbose,@options.debug)

    inputVars.each {|ivar| inVars[ivar] = nil}

    inVars
  end

  def getAndCheckOutputVariables
    outVars = {}
    oVars   = @options.outputVars.nil? ? Ecmwf2Icon::DefaultOutputVars : @options.outputVars

    Dbg.msg("Outout vars to be checked: #{oVars.join(" ")}",false,@options.debug)
    # * grid information must be provided in the configuration
    oVars.each {|outVar|
      Dbg.msg("Checking output variable '#{outVar}' for grid desc in configuration ... ",false, @options.debug)
      unless Ecmwf2Icon.hasGridinfo?(outVar,@config)
        warn "Cannot find grid information of variable '#{outVar}'!"
        warn "Abort #{__FILE__} !"
        exit 1
      else
        Dbg.msg("Output variable '#{outVar}' is valid and will be put on #{Ecmwf2Icon.grid(outVar,@config)} grid",@options.verbose,false)
        outVars[outVar] = nil
      end
    }
    outVars
  end

  # Split the input file into one file per variable
  def readVars
    threads = []
    @inVars.keys.each {|k|
      Dbg.msg("Reading variable '#{k}' ...",false,@options.debug)
      varfile = tfile
      threads << Thread.new(k) {|v| readVar(v,varfile) }
    }
    threads.each {|t| t.join}
  end

  # Every variable 'var' is read from the input and saved into another temporary file
  def readVar(var,varfile)
    operator = case var
               when Fixnum then 'selcode'
               when String then 'selname'
               else
                 warn "Wrong usage of variable identifier for '#{var}' (class #{var.class})!"
               end
    Cdo.send(operator,var,:in => @ifile, :out => varfile)
    @lock.synchronize {@inVars[var] = varfile}
  end

  # create a file for a given varriable on a certain grid
  def createVar(name); end

  # For which input variabes does exist a valid configuration? 
  # use names for netcdf input and codes for grib1
  def checkConfigOnInputVariables
    possibleOutputVars = []
    @inVars.each_key {|var|
      case var
      when Fixnum
        conf = @config.selectBy(:code => var)
      when String
        conf = @config.selectBy(:inputname => var)
      else
        warn "Wrong class for variable identifier of '#{var}'!"
        exit 1
      end
      if conf.size != 0
        outvar, invar, code, desc = conf.datasets(:outputname, :inputname, :code, :Description)[0]
        possibleOutputVars << outvar
        Dbg.msg("Found info for code:#{code} => dwd:#{outvar}, ecmwf:#{invar}, desc:#{desc})",
                @options.verbose,
                @options.debug)
      else
        Dbg.msg("Input variable with code '#{var}' has no config",
                @options.verbose,
                @options.debug)
      end
    }
    possibleOutputVars
  end

  # Main method for setting the preprocessing output
  def computeOutputVars
    Dbg.msg("Start computing output variables",false,@options.debug)
    Dbg.msg("Do only remapping of the required output varialbes if they are present in the input file",
            @options.verbose,@options.debug)

    ths = []
    @outVars.each_key {|ovar|
      Dbg.msg("Processing '#{ovar}' ...",@options.verbose,@options.debug)

      ths << Thread.new(ovar) {|ovar|
      # determine the grid of the variable
      grid     = Ecmwf2Icon.grid(ovar,@config).to_sym
      gridfile = @gridVars[GRID_VARIABLES[grid]]

      # remap the variable onto its icon grid provided in the configuration
      ivar              = @config.selectBy(:outputname => ovar).code[0].to_i
      copyfile, outfile = tfile, tfile

      # create netcdf version of the input
      Cdo.chainCall("setname,#{ovar} -copy",:in => @inVars[ivar],:out => copyfile,:options => "-f nc")

      # Perform conservative remapping
      Cdo.remapcon(gridfile,:in => copyfile,:out => outfile,:options => "-P #{@options.openmp}")

      @lock.synchronize { @outVars[ovar] = outfile }
      }
    }
    ths.each {|t| t.join}
  end

  def applyLSM
  end

  def consistencyCheck
  end

  def writeOutput(ofile='ifs2icon.grb')
    Dbg.msg("Compining files into output file '#{ofile}'",@options.verbose, @options.debug)
    Cdo.merge(:in => @outVars.values.join(" "),:out => ofile)
  end

  def remapVar
  end

  def iCodes
    @inVars.keys
  end
  def tfile
    t = Tempfile.new(self.class.to_s)
    @@_tempfiles << t
    t.path
  end
  private :readVars, :readVar, :remapVar, :tfile
end

# MEMO ========================================================================
#  grid description is done via the folling shortcuts:
#   |v|vertex|
#   |c|cell|
#   |e|edge|
#   |t|time|
#   |l|level|
# a variabel with lifes on the a (time,level,cell) grid wil get the description tlc
#
# ATMOSPHERE:
# Hydrostatic states: 
# TYPE t_hydro_atm_prog
#   REAL(wp), ALLOCATABLE ::  &
#   & pres_sfc(:,  :),  &!< surface pressure [Pa]        (nproma,     nblks_c)
#   &       vn(:,:,:),  &!< normal wind [m/s]            (nproma,nlev,nblks_e)
#   &     temp(:,:,:),  &!< temperature [K]              (nproma,nlev,nblks_c)
#   &    theta(:,:,:),  &!< potential temperature [K]    (nproma,nlev,nblks_c)
#   &   tracer(:,:,:,:)  !< tracer concentration [kg/kg] (nproma,nlev,nblks_c,ntracer)
# END TYPE t_hydro_atm_prog
# TYPE t_hydro_ocean_prog
#
#   REAL(wp), ALLOCATABLE ::    &
#     &  h(:,:) [ELEV]               ,& ! height of the free surface. Unit: [m] ! dimension:(nproma, nblks_c)
#     &  vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s] ! dimension: (nproma, n_zlev, nblks_e) ICON can compute vn out of u and v (primal_map_e2c)
#     &  tracer(:,:,:,:)          ! tracer concentration.  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce)
#                                 ! Ordering of tracers: 
#                                 !   1) pot_temp:= potential temperature, Unit: [deg C] (pt)
#                                 !   2) salinity:= salinity, Unit [psu] (sal)
#
# END TYPE t_hydro_ocean_prog
#
# OCEAN:
# MPIOM:
# CALL varlist_add(ocean_varlist,new_var('zo','sea_surface_height_above_sea_level','m',1,zo,'p','s'))
# CALL varlist_add(ocean_varlist,new_var('tho','sea_water_potential_temperature','C',2,tho,'p','c'))
# CALL varlist_add(ocean_varlist,new_var('uko','sea_water_x_velocity','m s-1',3,uko,'u','c',staggered=.true.))
# CALL varlist_add(ocean_varlist,new_var('vke','sea_water_y_velocity','m s-1',4,vke,'v','c',staggered=.true.))
# CALL varlist_add(ocean_varlist,new_var('sao','sea_water_salinity','psu',5,sao,'p','c'))
# =============================================================================

if __FILE__ == $0
  test = ENV['test'].nil? ? 'proc' : ENV['test']
  case test
  when 'proc'
    options  = PreProcOptions.parse(ARGV)
    pp options if options.debug
    exit
    Cdo.Debug = options.debug if options.verbose
    #=======================================================
    p    = Ifs2Icon.new(options)
    if options.check
      p.checkForPossibleOutput
    else
      p.run
    end
  when 'ext'
    conf = Ecmwf2Icon::DefaultConfig
    puts Ecmwf2Icon::DefaultConfig.to_string("tsv")
    #   puts Ecmwf2Icon::DefaultConfig.datacolumns
    puts "ecmwf's name of 'OMEGA' is '#{conf.selectBy(:outputname => 'OMEGA').inputname[0]}'"
  when 'cdo'
    include Cdo
    Cdo.Debug = ENV['debug'].nil? ? false : true
    ifile = '/home/ram/data/examples/T.jan.nc'
    puts Cdo.seltimestep(1,:in => Cdo.selname('T',:in =>ifile),:out => 'ofile.nc')
    pp Cdo.showcode(:in => ifile)
    pp Cdo.showcode(:in => 'ofile.nc')
    pp Cdo.showname(:in => "IFS.grb")
  end
end
