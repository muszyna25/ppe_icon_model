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

  # Only operators with documentation are accessible vie the build-in help.
  # Other have to be added manually
  @@undocumentedOperators = %w[geopotheight pressure_fl pressure_hl]

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
# ==============================================================================
# Option handling
class PreProcOptions
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    options                        = {}
    options[:interpolation_type]   = :simple
    options[:verbose]              = false
    options[:openmp]               = 4
    options[:model_type]           = 'hydrostatic'
    options[:threaded]             = false
    options[:persistent_tempfiles] = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: ifs2icon.rb [options]"

      opts.separator ""
      opts.separator "Specific options:"

      # Mandatory arguments
      #   global gridfile
      opts.on("-g", "--grid-file FILE",
              "Use FILE to read grid definition/resolution") do |file|
        options[:gridfile] = file
              end
      #   or preselected cell grid
      opts.on("--cell-grid-file FILE",
              "Use FILE for remapping onto the cell grid (precomputed weights)") do |file|
        options[:cellgridfile] = file
              end
      #   or precomputed weight files
      opts.on("--cell-weight-file FILE",
              "Use FILE for remapping onto the cell grid (precomputed weights)") do |file|
        options[:cellweightfile] = file
              end
      opts.on("--edge-weight-file FILE",
              "Use FILE for remapping onto the edge grid (precomputed weights)") do |file|
        options[:edgeweightfile] = file
              end

      opts.on("-i", "--input-file FILE",
              "Use FILE read input variables") do |file|
        options[:inputfile] = file
              end
      opts.on("-o", "--output-file FILE",
              "Use FILE write output variables") do |file|
        options[:outputfile] = file
              end
      opts.on("-O", "--output-variables x,y,z", Array, "'list' of output variables") do |list|
        options[:outputvars] = list
      end

      # Optional argument with keyword completion.
      opts.on("-m", "--model-type [TYPE]", [:hydrostatic, :nonhydrostatic],
              "Select model type (hydrostatic, nonhydrostatic) - default:hydrostatic") do |t|
        options[:model_type] = t
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

      opts.on("-V","--vct-file FILE","FILE should containt a vertable coordinates table for the required number of z levels") do |file|
        options[:vctfile] = file
      end
      opts.on("-O","--orography-file FILE","FILE should containt the required orography") do |file|
        options[:orofile] = file
      end

      # Optional argument with keyword completion.
      opts.on("-t", "--interpolation-type [TYPE]", [:simple, :bilinear, :bicubic, :horizontal_only],
              "Select interpolation type (simple,bilinear,bicubic,horizontal_only)") do |t|
        options[:interpolation_type] = t
      end

      # Set the path of the cdo binary to use
      opts.on("--cdo PATH_TO_CDO","Set the path of the CDO binary for computations (default:/usr/bin/cdo)") {|path|
        options[:cdo] = path
      }

      # Set OpenMP multithreadding for CDO
      opts.on("-P", "--openmp p", Numeric, "Set number of OpenMP threads to <p>") do |p|
        options[:openmp] = p
      end
      # Switch on multithreadded processing files/varialbles
      opts.on("-T","--multithreaded","Process files/variables/remapping multithreaded (Can produse high load!)") do |v|
        options[:threaded] = v
      end

      # create zlevel options for each given output variables
      opts.on("-z","--zlevels lev0,lev1,..,levN",
              Array,
              "'list' of target zlevels (not required for initial ICON data)") {|list|
        options[:zlevels] = list
      }
      opts.on("-Z","--zcoordinates FILE",
              String,
              "File with the 3D height field on full levels (nonhydrostatic case)") {|file|
        options[:fullzcoordinates] = file
      }
      opts.on("-Z","--full-zcoordinates FILE",
              String,
              "File with the 3D height field on full levels (nonhydrostatic case)") {|file|
        options[:fullzcoordinates] = file
      }
      opts.on("-Z","--half-zcoordinates FILE",
              String,
              "File with the 3D height field on half levels (nonhydrostatic case)") {|file|
        options[:halfzcoordinates] = file
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
        puts 'ifs2icon ' + [1,1].join('.')
        puts <<-END
        Performs horizontal and vertical interpolation off half and full level variables.
        Special handling for
        * pressure:
          compute differnce sto hydrostatic pressure using cdo stdatm operator
          perform vertical interpolation off these differences
          add hydrostatic pressure of the target vertical grid
        * vertical velocity
          input is OMEGA [Pa/s], output should be W [m/s]: OMEGA = -rho*g*W
          compute rho = P/(R*T) with P and T being hydrostatic (again cdo's stdatm)
          replace input vertical velocity with W=OMEGA/(-rho*g)
          END
        exit
      end

      # development switches: USE AT YOUR OWN RIST!
      opts.on("--persistent-tempfiles",
              "DEVELOPEMENT SWITCH: Create persistent intermediate files instead of cleaning up when finishing the script") {|v|
        options[:persistent_tempfiles] = v
      }

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
    end

    opts.parse!(args)

    # Consistency checks for options
    # mandatory ones
    begin
      mandatory = [:outputfile, :inputfile, :gridfile, :outputvars]
      missing = mandatory.select{ |param| options[param].nil? }
      if not missing.empty?
        warn "Missing options: #{missing.join(', ')}"
        puts opts
        exit
      end
      # test on file presence
      [:inputfile,:gridfile,:configfile].select {|param| not options[param].nil? }.each {|tag|
        unless File.exist?(options[tag])
          warn "Could not file #{tag.to_s} #{options[tag]}!"
          exit -1
        end
      }
    rescue OptionParser::InvalidOption, OptionParser::MissingArgument
      puts $!.to_s
      puts opts
      exit
    end

    checkModelOptions(options)

    options
  end
  def self.splitLevels(list)
    ret = {}
    list.slice_before {|v| v.to_i.to_s != v}.each {|sublist|
      ret[sublist[0]] = sublist[1..-1]
    }
    ret
  end
  def self.checkFile(filename,sinfo,relatedOption,model)
    unless filename
      warn "Please provide #{sinfo} file (option:#{relatedOption}) for #{model} model input"
      exit
    end
    unless File.readable?(filename) and File.file?(filename)
      warn "Unable to read from #{sinfo} file '#{filename}'"
      exit
    end
  end
  def self.checkModelOptions(options)
    pp options
    unless options[:interpolation_type].to_sym == :horizontal_only
      # for hydrostatic mode vct and orography are required
      if options[:model_type].to_s == "hydrostatic"
        checkFile(options[:vctfile], 'vct'      , '--vct-file'      , options[:model_type])
        checkFile(options[:orofile], 'orography', '--orography-file', options[:model_type])
        # check vct file: Should have C-style notation (start with 0) which is not
        # the case for icon vcts
        if options[:vctfile]
          uncommentedLines = File.open(options[:vctfile]).readlines.map(&:chomp).find_all {|l| l[0] != '#'}.map(&:strip)
          # check if first line starts with 1 or 0; if 1 increate all indicess by 1
          if uncommentedLines[0][0] == "1"
            newvct = options[:vctfile] + '_converted'
            Dbg.msg("VCT start with 1 but CDO requires C-Style indexing. Create converted vct file '#{newvct}'",options[:verbose],options[:debug])
            File.open(newvct,"w") {|f|
              uncommentedLines.each {|line| i,a,b = line.gsub(/ +/,' ').split; f << [i.to_i-1 ,a,b].join(' ') << "\n"}
            }
            options[:vctfile] = newvct
          end
        end
      end
      # target vertical coordinate is required for nonhydrostatic mode
      if options[:model_type].to_s == "nonhydrostatic"
        # check for presence of 3d vertical coordinate on full level if there is such an output variable
        levelKey = Ecmwf2Icon::DefaultColumnNames[:nlevels]
        unless (options[:config].send(levelKey).uniq & Ecmwf2Icon::DefaultNlevelIDs[:full]).empty?
          checkFile(options[:fullzcoordinates], 'target vertical coordinate', '--full-zcoordinates', options[:model_type]) 
        end
        unless (options[:config].send(levelKey).uniq & Ecmwf2Icon::DefaultNlevelIDs[:half]).empty?
          checkFile(options[:halfzcoordinates], 'target vertical coordinate', '--half-zcoordinates', options[:model_type]) 
        end
      end
    end
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

  DefaultColumnNames = {
      :outName                 => :outputname,
      :inName                  => :inputname,
      :code                    => :code,
      :grid                    => :grid,
      :levType                 => :typeOfLayer,
      :nlevels                 => :nlevel,
      :info                    => :description,
      :horizontalInterpolation => :I
  }

  DefaultNlevelIDs = {
    :full => ['full','nlev'],
    :half => ['half','nlev+1'],
    :surf => ['surf','1']
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
  C_R          = 287.05
  SCALEHEIGHT  = 10000.0
  C_EARTH_GRAV = 9.80665
  # function for later computation of hydrostatic atmosphere pressure
  PRES_EXPR = lambda {|height| "101325.0*exp((-1)*(1.602769777072154)*log((exp(#{height}/#{SCALEHEIGHT})*213.15+75.0)/288.15))"}
  TEMP_EXPR = lambda {|height| "213.0+75.0*exp(-#{height}/#{SCALEHEIGHT})"}
  RHO_EXPR  = lambda {|pressure,temperature| "#{pressure}/(#{C_R}*#{temperature})"}
  W_EXPR    = lambda {|omega,rho| "#{omega}/(-#{rho}*#{C_EARTH_GRAV})"}

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

    @invars, @outvars, @grids, @weights = {}, {}, {}, {}

    @ifile, @ofile = @options[:inputfile], @options[:outputfile]

    @itype = '.grb' == File.extname(@ifile) ? :code : :name

    # Create internal configuration state
    if @options[:config].nil?
      warn "Configuration is not provided. Internal default is used instead."
      @config = Ecmwf2Icon::DefaultConfig
    else
      @config = @options[:config]
    end

    Ecmwf2Icon.checkConfig(@config)

    Ecmwf2Icon.displayConfig(@config,@options[:verbose], @options[:debug])

    prepareGridAndWeights


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
        conf = @config.selectBy(:inputname => var.upcase)
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
    splitInput
    #readVars

    Dbg.msg("Start computing output variables",false,@options[:debug])
    Dbg.msg("Do only remapping of the required output varialbes if they are present in the input file",
            @options[:verbose],@options[:debug])

    horizontalInterpolation

    if @options[:interpolation_type].to_s == "horizontal_only"
      # merge all horizontal interpolations together
      Dbg.msg("Creating output with horizontal interpolation only!",true,true)
      Dbg.msg("Combining files into output file '#{@ofile}'",@options[:verbose], @options[:debug])
      Cdo.merge(:in => @outvars.values.join(" "),:out => @ofile)
      exit
    else
      verticalInterpolation
    end

    applyLSM

    consistencyCheck

    writeOutput
  end

  # Create files for later usage as grid description for the remapping
  def prepareGridAndWeights
    GRID_VARIABLES.each {|k,v|
      # only generate weights if they are needed
      next unless @config.grid.uniq.include?(k.to_s)

      if k == :cell and not @options[:cellgridfile].nil?
        Dbg.msg("Use #{k} grid from #{@options[:cellgridfile]}",false, @options[:debug])
        @grids[v] =  @options["#{k}gridfile".to_sym]
      else
        Dbg.msg("Selecting #{k} grid from #{@options[:gridfile]}",false, @options[:debug])
        varfile = tfile
        Cdo.selname(v,:in => @options[:gridfile],:out => varfile)
        @grids[v] = varfile
        Dbg.msg("#{k} grid is now in #{varfile}",false, @options[:debug])
      end

      unless @options[:cellweightfile].nil?
        Dbg.msg("Using #{k} weights from precomputed file: #{@options[:cellweightfile]}",false, @options[:debug])
        @weights[v] = @options["#{k}weightfile".to_sym]
      else
        weightsfile = @options[:debug] ? "#{k.to_s}_weight.nc" : tfile
        Dbg.msg("Compute #{k} weights",false, @options[:debug])
        Cdo.gencon(varfile,:in => @ifile, :out => weightsfile,:options => "-P #{@options[:openmp]}")
        @weights[v] = weightsfile
      end
    }
  end

  def getInputVariables
    invars = {}
    # Use names for netcdf files and codes for grib input files
    operator  = @itype == :code ? :showcode : :showname
    inputVars = Cdo.send(operator,:in => @ifile)

    # Use numbers as identifiers for grib input
    inputVars.map!(& :to_i) if operator == :showcode

    Dbg.msg("Found these variables in the input file = #{inputVars.join(",")} in input file",@options[:verbose],@options[:debug])

    inputVars.each {|ivar| invars[ivar] = nil}

    invars
  end

  def fullLevelInputVars
    inputVarsByLevelType(:full)
  end
  def halfLevelInputVars
    inputVarsByLevelType(:half)
  end
  def surfLevelInputVars
    inputVarsByLevelType(:surf)
  end

  def inputVarsByLevelType(ltype)
    vnames = []
    tag = case ltype
    when :full then 'nlev'
    when :half then 'nlev+1'
    when :surf then '1'
    end
    @invars.each_key {|k|
      column = case k
      when Fixnum then :code
      when String then :inputname
      else 
        warn "Wrong Input variable information in variable 'invars'"
        exit
      end
      unless ( names = @config.selectBy(column => k,:nlevel => "'#{tag}'").outputname; names.empty?)
        vnames << names[0]
      end
    }
    vnames
  end

  def levelTypeOfVar(varname)
    [:full, :half, :surf].each {|ltype|
      return ltype if inputVarsByLevelType(ltype).include?(varname)
    }
    warn "No level Type found for variable '#{varname}'"
    return false
  end


  # Split the input file into one file per variable
  def readVars
    if @options[:threaded]
      threads = []
      @invars.keys.each {|k|
        Dbg.msg("reading variable '#{k}' ...",false,@options[:debug])
        varfile = tfile
        threads << Thread.new(k) {|var| 
          readVar(var,varfile)
          @lock.synchronize {@invars[var] = varfile}
        }
      }
      threads.each {|t| t.join}
    else
      @invars.keys.each {|var|
        Dbg.msg("reading variable '#{var}' ...",false,@options[:debug])
        varfile = tfile
        readVar(var,varfile)
        @invars[var] = varfile
      }
    end
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
  end

  def splitInput(tag=Time.new.strftime("%Y%m%d-%H%M%S")+"_")
    operator  = @itype == :code ? :splitcode : :splitname
    Cdo.send(operator,:in => @ifile,:out => tag)
    splitfiles = Dir.glob("#{tag}*.*")

    cmp = @itype == :code \
      ? lambda {|tagFromFile,varID| tagFromFile.to_i == varID} \
      : lambda {|tagFromFile,varID| tagFromFile == varID}
  
    ids = splitfiles.map {|f| f.split('_')[-1].split('.')[0]}
    @invars.each_key {|var|
      @invars[var] = splitfiles.find_all {|f| cmp.call(f.split('_')[-1].split('.')[0],var) }.first
    }

    pp @invars
  end

  def horizontalInterpolation4Var(outvar)
    Dbg.msg("Processing '#{outvar}' ...",@options[:verbose],@options[:debug])
    # determine the grid of the variable
    grid       = Ecmwf2Icon.grid(outvar,@config).to_sym
    gridfile   = @grids[GRID_VARIABLES[grid]]
    weightfile = @weights[GRID_VARIABLES[grid]]

    # remap the variable onto its icon grid provided in the configuration
    meth = @itype == :code ? :code : :inputname
    ivar = @itype == :code ? @config.selectBy(DefaultColumnNames[:outName] => outvar).code[0].to_i \
                           : @config.selectBy(DefaultColumnNames[:outName] => outvar).inputname[0].downcase
    copyfile, outfile = tfile, tfile

    # Perform conservative remapping with pregenerated weights and rename to target variable names
    Cdo.chainCall("-remap,#{gridfile},#{weightfile} -setname,#{outvar} -copy",
                  :in => @invars[ivar],
                  :out => outfile,
                  :options => "-f nc2")

    outfile
  end
  def horizontalInterpolation
    if @options[:threaded]
      ths = []
      @outvars.each_key {|outvar|
        ths << Thread.new(outvar) {|ovar| 
          outfile = horizontalInterpolation4Var(ovar)
          @lock.synchronize { @outvars[ovar] = outfile }
        }
      }
      ths.each {|t| t.join}
    else
      @outvars.each_key {|ovar|
        outfile = horizontalInterpolation4Var(ovar)
        @outvars[ovar] = outfile
      }
    end
  end

  def verticalInterpolation
    # merge all horizontal interpolations together
    intermediateFile, hybridlayerfile, reallayerfile = tfile, tfile, tfile
    Dbg.msg("Combining files into output file '#{intermediateFile}'",false, @options[:debug])
    Dbg.msg("Output variables before vertical interpolation: " + @outvars.keys.join(" "),false, @options[:debug])
    Cdo.merge(:in => @outvars.values.join(" "),:out => intermediateFile)

    case @options[:model_type]
    when 'hydrostatic'
      # perform vertical interpolation wrt. original surface pressure and orography
      Cdo.remapeta(@options[:vctfile],@options[:orofile],:in => intermediateFile,:out => hybridlayerfile)
      unless @options[:zlevels].nil?
        # perform hybrid2realLevel conversion
        Cdo.ml2hl(@options[:zlevels].reverse.join(','),:in => hybridlayerfile, :out => reallayerfile)
        @_preout = reallayerfile
      else
        @_preout = hybridlayerfile
      end
    when 'nonhydrostatic'
      # create 3d height from IFS intput data with 'cdo geopotheight'
      intermediateHeight = tfile
      Cdo.geopotheight(:in => intermediateFile,:out => intermediateHeight)

      # Create output with vertical height axis of intermediate IFS (debug only)
      Cdo.copy(:in => intermediateHeight,:out => 'intermediateHeights.nc') if @options[:debug]
      Cdo.copy(:in => intermediateFile,  :out => 'intermediateIFS.nc')     if @options[:debug]

      iconPressure     = verticalPressureHandling(intermediateFile,intermediateHeight)
      intermediateFile = verticalVelocityHandling(intermediateFile,intermediateHeight)

      # =======================================================================
      # HANDLING OF OTHER VARIABLES
      # * perform 3d interpolation of any other 3d variables with intlevel3d
      remainingVariablesFile, @_preout = tfile,tfile
      # split into full and half level variables
      fullLevelVarnames = fullLevelInputVars
      halfLevelVarnames = halfLevelInputVars
      surfLevelVarnames = surfLevelInputVars
      pp({:full => fullLevelVarnames, :half => halfLevelVarnames,:surf => surfLevelVarnames}) if @options[:debug]

      iconOutFiles = []
      {
        :full => fullLevelVarnames,
        :half => halfLevelVarnames,
        :surf => surfLevelVarnames
      }.each {|ltype,lvarnames|
        unless lvarnames.empty?
          selfile = tfile
          Cdo.selname(lvarnames.join(','),:in => intermediateFile,:out => selfile)
          targetZFile = case ltype 
                        when :full then @options[:fullzcoordinates]
                        when :half then @options[:halfzcoordinates]
                        end
          unless ltype == :surf
            iconLevelVarsFile = tfile
            Cdo.intlevelx3d(intermediateHeight, :in => [selfile,targetZFile].join(' '), :out => iconLevelVarsFile)
            iconOutFiles << iconLevelVarsFile
          else
            iconOutFiles << selfile
          end
        else
          puts "No #{ltype.to_s} type variables found." if @options[:debug]
        end
      }


      # =======================================================================
      # MERGE BACK INTO A SINGLE OUTPUT FILE
      Cdo.merge(:in => [iconOutFiles,iconPressure].flatten.join(' '),:out => @_preout)
    end
  end

  def verticalPressureHandling(intermediateFile,intermediateHeight)
    # some variable name definitions. They might change in the future
    iFShydpressVarname         = 'pres'
    iCONhydpressVarname        = 'pres'
    iCONverticalCoordinateName = 'ZF3'

    # =======================================================================
    # PRESSURE HANDLING:
    # * create 3d pressure from intermediate IFS interpolation (horiz)
    intermediatePressure = tfile
    Cdo.chainCall("-setname,#{iFShydpressVarname} -pressure_fl",:in => intermediateFile,:out => intermediatePressure)

    # debug output containing the intermediate 3d IFS pressure coordinate
    Cdo.copy(:in => intermediatePressure,:out => 'intermediateIFSPressure.nc') if @options[:debug]

    if false #version for patched geopotheight operator
      Cdo.selname('pres',:in => intermediateHeight, :out => intermediatePressure)
      tmp = tfile
      Cdo.delname('pres',:in => intermediateHeight, :out => tmp)
      Cdo.copy(:in => tmp,:out => intermediateHeight)
    end

    # * substract hydrostatic pressure on intermediate IFS heights from the
    #   real intermediate IFS pressure levels
    #   - compute hydrostatic pressure on 3d height field, same formular like stdatm operator in CDO
    #     * Use expr operator like
    #        cdo expr 'press=101325.0*exp((-1)*(1.602769777072154)*log((exp(height/10000.0)*213.15+75.0)/288.15))'
    hydrostaticPresOnIntermediateHeights, diff, diffOnIconVertGrid = tfile,tfile,tfile
    Cdo.expr("'#{iFShydpressVarname}=#{PRES_EXPR['geopotheight']}'",
             :in => intermediateHeight, :out => hydrostaticPresOnIntermediateHeights)

    #     * substract the hydrostatic pressure from the IFS pressures
    Cdo.sub(:in => [intermediatePressure,hydrostaticPresOnIntermediateHeights].join(' '), :out => diff)
    #
    # * 3d linerar interpolation if the pressure differences onto the ICON vertical height levels
    Cdo.intlevelx3d(intermediateHeight,:in => [diff,@options[:fullzcoordinates]].join(' '), :out => diffOnIconVertGrid)
    #
    # * add hydrostatic pressure values for the ICON heights
    #   - remove target z coordinate and reset variable name for later operations (add)
    tmp = tfile
    Cdo.delname(iCONverticalCoordinateName, :in => diffOnIconVertGrid, :out => tmp)
    Cdo.setname(iCONhydpressVarname       , :in => tmp               , :out => diffOnIconVertGrid) #TODO not necessary
    # prepare hydrostatic pressure in vertical ICON grid
    hydrostaticPressOnIconHeights, iconPressure = tfile, tfile
    Cdo.expr("'#{iCONhydpressVarname}=#{PRES_EXPR[iCONverticalCoordinateName]}'",
             :in => @options[:fullzcoordinates], :out => hydrostaticPressOnIconHeights)
    Cdo.add(:in => [diffOnIconVertGrid,hydrostaticPressOnIconHeights].join(' '),:out => iconPressure)

    # Interpolated pressure for later debugging
    Cdo.copy(:in => iconPressure,:out => 'iCONPressure.nc') if @options[:debug]

    iconPressure
  end

  # Compute W[m/s] out of OMEGA [Pa/s] according to
  # * OMEGA = -rho*g*W
  # * rho   = P/(R*T) out of P/rho = R*T
  # where P and T are hydrostatic pressure and temperatur. After this
  # computation the vertical interpolation onto the target is done so that
  # there is only on interpolation step
  #
  # verticalVelocityHandling() changes the name of the intermediate File
  # (pre-vertical-interpolation-state)
  def verticalVelocityHandling(intermediateFile, intermediateHeight)
    # select variable from intermediate IFS intput AND REMOVE it from the origin
    verticalVelocityName = @config.selectBy(:code => 135).outputname[0]
    intermediateW, intermediateWWithNewUnit, tmp = tfile, tfile, tfile
    Cdo.selname(verticalVelocityName,:in => intermediateFile, :out => intermediateW)

    # create pressure and temperature for appropriate vertical coordinate
    presFile, tempFile, densityFile = tfile , tfile, tfile
    presName, tempName, densityName = 'pres', 't'  , 'rho'
    Cdo.expr("'#{presName}=#{PRES_EXPR['geopotheight']}'", :in => intermediateHeight, :out => presFile)
    Cdo.expr("'#{tempName}=#{TEMP_EXPR['geopotheight']}'", :in => intermediateHeight, :out => tempFile)
    # compute the hydrostatic density
    #Cdo.chainCall("setname,#{densityName} -mulc,#{C_R} -div",in: [presFile,tempFile].join(' '), out: densityFile)
    Cdo.chainCall("setname,#{densityName} -divc,#{C_R} -div",:in => [presFile,tempFile].join(' '), :out => densityFile)
    Cdo.merge(:in => [intermediateW,densityFile].join(' '),:out => tmp)
    # compute vertical velocity in new units
    Cdo.expr("'#{verticalVelocityName}=#{W_EXPR[verticalVelocityName,densityName]}'",:in => tmp, :out => intermediateWWithNewUnit)

    Cdo.copy(:in => intermediateWWithNewUnit, :out => 'intermediateVerticalVelocity.nc') if @options[:debug]


    # add the new vertical velocity to the intermediat file back again
    newIntermediateFile = tfile
    Cdo.replace(:in => [intermediateFile,intermediateWWithNewUnit].join(' '), :out => newIntermediateFile)

    newIntermediateFile
  end
  def writeOutput(ofile=@ofile)
    Dbg.msg("Create output file '#{ofile}'",@options[:verbose],@options[:debug])
    if @_preout
      Cdo.copy(:in => @_preout,:out => ofile)
    else
      Cdo.merge(:in => @outvars.values.join(' '), :out => ofile)
    end
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
    unless @options[:persistent_tempfiles]
      t = Tempfile.new(self.class.to_s)
      @@_tempfiles << t
      t.path
    else
      t = "_"+rand(10000000).to_s
      @@_tempfiles << t
      t
    end
  end
  private :readVars, :readVar, :remapVar, :tfile
end

if __FILE__ == $0
  test = ENV['test'].nil? ? 'proc' : ENV['test']
  case test
  when 'proc'
    #=======================================================
    # test on using Ruby 1.9.x
    if RUBY_VERSION < '1.9.0'
      warn "Please use Ruby 1.9.x! Your version is #{RUBY_VERSION}"
      exit -1
    end

    #=======================================================
    # options handling
    options  = PreProcOptions.parse(ARGV)
    Cdo.setCdo(options[:cdo]) if options[:cdo]
    pp options if options[:debug]
    Cdo.Debug = options[:debug] if options[:verbose]
    #=======================================================
    p    = Ifs2Icon.new(options)
    if options[:check]
      p.checkForPossibleOutput
    else
      p.run
    end
  when 'test'
    options  = PreProcOptions.parse(ARGV)
    Cdo.setCdo(options[:cdo]) if options[:cdo]
    pp options if options[:debug]
    Cdo.Debug = options[:debug] if options[:verbose]
    #=======================================================
    p    = Ifs2Icon.new(options)
    p.verticalVelocityHandling('','')
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
#
# vim:foldmethod=syntax
