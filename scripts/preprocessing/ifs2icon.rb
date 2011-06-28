#!/usr/bin/env ruby
# Due to the fact, that it will make heavy use of CDO, there is an CDO issue
# for that topic: https://code.zmaw.de/issues/612
#
# Main ICON issue is https://code.zmaw.de/issues/735
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
      warn "Operators could not get listed by running the CDO binary (#{@@CDO})"
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
    firstOp = chain[0...[chain.index(','),chain.index(' ')].min]
    firstOp = firstOp[1..-1] if firstOp[0] == '-'
    if /(info|show|griddes)/.match(firstOp)
      Cdo.run(" #{chain} #{io[:in]} ",$stdout)
    else
      opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
      Cdo.run(" #{chain}#{opts} #{io[:in]} ",io[:out],io[:options])
    end
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
end
# ==============================================================================
# Option handling
class PreProcOptions
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    options                    = {}
    options[:interpolation_type] = :simple
    options[:verbose]            = false
    options[:openmp]             = 4
    options[:model_type] = 'hydrostatic'

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: ifs2icon.rb [options]"

      opts.separator ""
      opts.separator "Specific options:"

      # Mandatory arguments
      opts.on("-g", "--grid-file FILE",
              "Use FILE read grid definition/resolution") do |file|
        options[:gridfile] = file
              end
      opts.on("-i", "--input-file FILE",
              "Use FILE read input variables") do |file|
        options[:inputfile] = file
              end
      opts.on("-o", "--output-file FILE",
              "Use FILE write output variables") do |file|
        options[:outputfile] = file
              end

      # optional configuration file using csv like format
      opts.on("--txt-config FILE", "Use FILE as configuration file", "Format should be a separated columns text file (e.g. csv)") do |file|
        options[:config] = Ecmwf2Icon::readConfig(file,'file')
      end
      # optional configuration file using JSON like format
      jsonStringExample=<<-EOF
                  [
                    ["description"       , "outputname", "inputname", "code", "grid", "typeOfLayer", "nlevel", "GP", "I", "Notes"],
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
      opts.on("--json-config FILE", "Use FILE as JSON configuration file","Format example:\n"+jsonStringExample) do |file|
        options[:config] = Ecmwf2Icon::readConfig(file,'file','json')
      end
      opts.on("-c", "--global-config", "<filename>", String, :REQUIRED) do |file|
        options[:configfile] = file
        # instead of using all these options, a config file can be specified
        load file #if Module.constants.include?("PreProcOpts")
        PreProcOpts.constants.each {|const|
          val = PreProcOpts.const_get(const)
          case const.to_s
          when 'CONFIG' then
            val = Ecmwf2Icon::readConfig(val,'string','json')
          when 'ZLEVELS' then
            val = val.split(',')
          when 'OUTPUTVARS' then
            val = val.split(',')
          else
            puts "Parsing value of #{const.to_s}" if options[:debug]
          end
          const = const.to_s.downcase.to_sym
          options[const] = val
        }
      end

      opts.on("-V","--vct-file FILE","FILE should containt a vertable coordinates table for the required number of z levels") do |file|
        options[:vctfile] = file
      end
      opts.on("-O","--orography-file FILE","FILE should containt the required orography") do |file|
        options[:orofile] = file
      end

      # Optional argument with keyword completion.
      opts.on("-t", "--interpolation-type [TYPE]", [:simple, :bilinear, :bicubic],
              "Select interpolation type (simple, bilinear, bicubic)") do |t|
        options[:interpolation_type] = t
      end
      # Optional argument with keyword completion.
      opts.on("-m", "--model-type [TYPE]", [:hydrostatic, :nonhydrostatic],
              "Select model type (hydrostatic, nonhydrostatic)") do |t|
        options[:model_type] = t
      end
      #
      # List of output variables
      opts.on("-O", "--output-variables x,y,z", Array, "'list' of output variables") do |list|
        options[:outputvars] = list
      end

      # Set OpenMP multithreadding for CDO
      opts.on("-P", "--openmp p", Numeric, "Set number of OpenMP threads to <p>") do |p|
        options[:openmp] = p
      end

      # Boolean switches
      opts.on("-s", "--strict", "Process only variables which have valid entry in the configuration") do |v|
        options[:strict] = v
      end
      opts.on("-C", "--check", "Check input file and configuration on which output variables are available.") do |v|
        options[:check] = v
      end
      opts.on("-v", "--verbose", "Run verbosely") do |v|
        options[:verbose] = v
      end
      opts.on("-D", "--debug", "Run in debug mode") do |v|
        options[:debug] = v
      end

      # create zlevel options for each given output variables
      opts.on("-z","--zlevels lev0,lev1,..,levN",
              Array,
              "'list' of target zlevels (not required for initial ICON data)") {|list|
        options[:zlevels] = list
      }
      options[:zcoordinates] = {}
      opts.on("-Z","--zcoordinates FILE",
              String,
              "File with the 3D height field (nonhydrostatic case)") {|file|
        options[:zcoordinates] = file
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
      mandatory = [:outputfile, :inputfile, :gridfile, :outputvars]
      missing = mandatory.select{ |param| options[param].nil? }
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
  def self.splitLevels(list)
    ret = {}
    list.slice_before {|v| v.to_i.to_s != v}.each {|sublist|
      ret[sublist[0]] = sublist[1..-1]
    }
    ret
  end
end
# ==============================================================================
# variable substitution between icon and Eecmwf
module Ecmwf2Icon
  _defaultConfig =[
    ['description'                        , 'outputname', 'inputname', 'code', 'grid', 'typeOfLayer', 'nlevel', 'GP', 'I', 'Notes']                                        ,
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
    ["geopotential "                      , "GEOSP"     , "GEOSP"    , "129" , "cell", "surface"    , "1"     , "1" , "Q", ""]     ,
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
      :info                    => :description,
      :horizontalInterpolation => :I
    }

  def defaultConfig
    DefaultConfig
  end

  def Ecmwf2Icon.readConfig(configData,dataType,dataFormat='txt')
    case dataType
    when 'file'
      case dataFormat
      when 'txt'
        config = ExtCsv.new("file","txt",configData)
        unless Ecmwf2Icon.check(@config)
          warn "Cannot read in file #{configFile}! Please check its structure"
          exit
        end
      when 'json'
        config = ExtCsv.new('array','plain',JSON.parse(File.open(configData,"r").read).transpose)
      else
        warn 'Unusable input type for configuration data'
        exit
      end
    when 'string'
      case dataFormat
      when 'json'
        config = ExtCsv.new('array','plain',JSON.parse(configData).transpose)
      end
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

  def number_of_processors
    if RUBY_PLATFORM =~ /linux/
      return `cat /proc/cpuinfo | grep processor | wc -l`.to_i
    elsif RUBY_PLATFORM =~ /darwin/
      return `sysctl -n hw.logicalcpu`.to_i
    elsif RUBY_PLATFORM =~ /(win32|mingw|cygwin)/
      # this works for windows 2000 or greater
      require 'win32ole'
      wmi = WIN32OLE.connect("winmgmts://")
      wmi.ExecQuery("select * from Win32_ComputerSystem").each do |system|
        begin
          processors = system.NumberOfLogicalProcessors
        rescue
          processors = 0
        end
        return [system.NumberOfProcessors, processors].max
      end
      elseif RUBY_PLATFORM =~ /java/
        return Runtime.getRuntime().availableProcessors().to_i
    end
    raise "can't determine 'number_of_processors' for '#{RUBY_PLATFORM}'"
  end
end

# A preprocessor should basically constist fof the following parts
# * knowledge about required in and output variables
# * rules for computing the output from the input variables
# * abtraction from the IO level
# The folloing steps have to be done
# * read input variables
# * compute the output variables
# * map output variables onto the required grid and vertical coordinates
# * (TODO) apply land-see-mask
# * (TODO) to concistency checks
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
  @@_preout    = ''

  attr_accessor :invars, :outvars, :config, :ifile

  def initialize(options)

    @lock = Mutex.new

    @options = options
    @options.freeze

    @tempfiles = []

    @invars, @outvars, @gridvars = {}, {}, {}

    # Create internal configuration state
    if @options[:config].nil?
      warn "Configuration is not provided. Internal default is used instead."
      @config = Ecmwf2Icon::DefaultConfig
    else
      @config = @options[:config]
    end

    Ecmwf2Icon.checkConfig(@config)

    Ecmwf2Icon.displayConfig(@config,@options[:verbose], @options[:debug])

    prepareGridfiles

    @ifile, @ofile = @options[:inputfile], @options[:outputfile]

    # Read input variables form File
    @invars = getInputVariables

    # Check which input/output variabels have a valid configuration
    checkConfigOnIOVariables

    @options[:outputvars].each {|ov| @outvars[ov] = nil}
  end

  # Check which input variable has a valid entry in the configuration and
  # ensure that the required output variables have valid grid info
  # After this check the list of output variables supplied by the user can be used
  def checkConfigOnIOVariables
    possibleOutputVars = []
    @invars.each_key {|var|
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
        outvar, invar, code, desc = conf.datasets(*DefaultColumnNames.values_at(:outName, :inName, :code, :info))[0]
        possibleOutputVars << outvar
        Dbg.msg("Found info for code:#{code} => dwd:#{outvar}, ecmwf:#{invar}, desc:#{desc})",
                @options[:verbose],
                @options[:debug])
      else
        Dbg.msg("Input variable with code '#{var}' has no config",
                @options[:verbose],
                @options[:debug])
      end
    }
    #check if all output variables have a config entry
    if @options[:outputvars] - possibleOutputVars != []
      warn "Outout variables #{(@options[:outputvars] - possibleOutputVars).join(',')} cannot be computed because of missing configuration"
      exit -1
    end
    #check if all output variables have grid information in config
    @options[:outputvars].each {|outVar|
      Dbg.msg("Checking output variable '#{outVar}' for grid desc in configuration ... ",false, @options[:debug])
      unless Ecmwf2Icon.hasGridinfo?(outVar,@config)
        warn "Cannot find grid information of variable '#{outVar}'!"
        warn "Abort #{__FILE__} !"
        exit 1
      else
        Dbg.msg("Output variable '#{outVar}' is valid and will be put on #{Ecmwf2Icon.grid(outVar,@config)} grid",@options[:verbose],false)
      end
    }
  end

  def run
    readVars

    computeOutputVars

    applyLSM

    consistencyCheck

    writeOutput
  end

  # Create files for later usage as grid description for the remapping
  def prepareGridfiles
    GRID_VARIABLES.each {|k,v|
      Dbg.msg("Selecting #{k} grid from #{@options[:gridfile]}",false, @options[:debug])
      varfile = tfile
      Cdo.selname(v,:in => @options[:gridfile],:out => varfile)
      @gridvars[v] = varfile
      Dbg.msg("#{k} grid is now in #{varfile}",false, @options[:debug])
      Dbg.msg(IO.popen("ls -crtlh #{varfile}").read,false,@options[:debug])
    }
  end

  def getInputVariables
    invars = {}
    # Use names for netcdf files and codes for grib input files
    operator  = /grb/.match(File.extname(@ifile)) ? :showcode : :showname
    inputVars = Cdo.send(operator,:in => @ifile)
    # Use numbers as identifiers for grib input
    inputVars.map!(& :to_i) if operator == :showcode

    Dbg.msg("Found these variables in the input file = #{inputVars.join(",")} in input file",@options[:verbose],@options[:debug])

    inputVars.each {|ivar| invars[ivar] = nil}

    invars
  end

  # Split the input file into one file per variable
  def readVars
    threads = []
    @invars.keys.each {|k|
      Dbg.msg("Reading variable '#{k}' ...",false,@options[:debug])
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
    @lock.synchronize {@invars[var] = varfile}
  end

  # Main method for setting the preprocessing output
  def computeOutputVars
    Dbg.msg("Start computing output variables",false,@options[:debug])
    Dbg.msg("Do only remapping of the required output varialbes if they are present in the input file",
            @options[:verbose],@options[:debug])

    horizontalInterpolation

    verticalInterpolation
  end

  def horizontalInterpolation
    ths = []
    @outvars.each_key {|ovar|
      Dbg.msg("Processing '#{ovar}' ...",@options[:verbose],@options[:debug])

      ths << Thread.new(ovar) {|ovar|
        # determine the grid of the variable
        grid     = Ecmwf2Icon.grid(ovar,@config).to_sym
        gridfile = @gridvars[GRID_VARIABLES[grid]]

        # remap the variable onto its icon grid provided in the configuration
        ivar              = @config.selectBy(DefaultColumnNames[:outName] => ovar).code[0].to_i
        copyfile, outfile = tfile, tfile

        # create netcdf version of the input
        Cdo.chainCall("-setname,#{ovar} -copy",:in => @invars[ivar],:out => copyfile,:options => "-f nc")

        # Perform conservative remapping
        Cdo.remapcon(gridfile,:in => copyfile,:out => outfile,:options => "-P #{@options[:openmp]}")

        @lock.synchronize { @outvars[ovar] = outfile }
      }
    }
    ths.each {|t| t.join}
  end

  def verticalInterpolation
    # merge all horizontal interpolations together
    intermediateFile, hybridlayerfile, reallayerfile = tfile, tfile, tfile
    Dbg.msg("Compining files into output file '#{intermediateFile}'",@options[:verbose], @options[:debug])
    Dbg.msg(@outvars.keys.join(" "),@options[:verbose], @options[:debug])
    Cdo.merge(:in => @outvars.values.join(" "),:out => intermediateFile)
    if @options[:model_type] == 'hydrostatic'
      # perform vertical interpolation wrt. original surface pressure and orography
      Cdo.remapeta(@options[:vctfile],@options[:orofile],:in => intermediateFile,:out => hybridlayerfile)
      unless @options[:zlevels].nil?
        # perform hybrid2realLevel conversion
        Cdo.ml2hl(@options[:zlevels].reverse.join(','),:in => hybridlayerfile, :out => reallayerfile)
        @_preout = reallayerfile
      else
        @_preout = hybridlayerfile
      end
    else
      #perform interpolation of 3D height field
      warn "Vertical interpolation onto a 3D vertical coordinate is not implemented, yet!"
      warn "Abort #{__FILE__}!"
      exit
    end
  end

  def writeOutput
    Cdo.copy(:in => @_preout,:out => @ofile)
  end

  # create a file for a given varriable on a certain grid
  def createVar(name); end
  def remapVar; end
  def applyLSM; end
  def consistencyCheck; end

  def iCodes
    @invars.keys
  end
  def tfile
    t = Tempfile.new(self.class.to_s)
    @@_tempfiles << t
    t.path
  end
  private :readVars, :readVar, :remapVar, :tfile
end

if __FILE__ == $0
  test = ENV['test'].nil? ? 'proc' : ENV['test']
  case test
  when 'proc'
    options  = PreProcOptions.parse(ARGV)
    pp options if options[:debug]
    Cdo.Debug = options[:debug] if options[:verbose]
    #=======================================================
    p    = Ifs2Icon.new(options)
    if options[:check]
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

# =============================================================================
# Open Points:
# * verification to hydrostatic model
# * orography on edges
# * 3D vertical coordinate
