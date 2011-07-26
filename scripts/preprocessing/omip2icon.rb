#!/usr/bin/env ruby
$LOAD_PATH << '.'

require 'tempfile'
require 'thread'
require 'pp'
require 'extcsv'
require 'optparse'
require 'json'

require 'ifs2icon'

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
class PreProc4Icon

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

    # Create internal configuration state
    if @options[:config].nil?
      warn "Configuration is not provided. Abort"
      exit
    else
      @config = @options[:config]
    end

    #TODO @config.check
    #TODO @config.display if @options[:verbose] or @options[:debug]

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

      Dbg.msg("Selecting #{k} grid from #{@options[:gridfile]}",false, @options[:debug])
      varfile = tfile
      Cdo.selname(v,:in => @options[:gridfile],:out => varfile)
      @grids[v] = varfile
      Dbg.msg("#{k} grid is now in #{varfile}",false, @options[:debug])

      weightsfile = tfile
      Cdo.gencon(varfile,:in => @ifile, :out => weightsfile,:options => "-P #{@options[:openmp]}")
      @weights[v] = weightsfile
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
      @invars.keys.each {|k|
        Dbg.msg("reading variable '#{k}' ...",false,@options[:debug])
        varfile = tfile
        readVar(k,varfile)
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

  def horizontalInterpolation
    if @options[:threaded]
      ths = []
      @outvars.each_key {|outvar|
        Dbg.msg("Processing '#{outvar}' ...",@options[:verbose],@options[:debug])

        ths << Thread.new(outvar) {|ovar|
          # determine the grid of the variable
          grid       = @config.selectBy(:outName => ovar).grid[0].to_sym
          gridfile   = @grids[GRID_VARIABLES[grid]]
          weightfile = @weights[GRID_VARIABLES[grid]]

          # remap the variable onto its icon grid provided in the configuration
          ivar              = @config.selectBy(:outName => ovar).inName[0]
          copyfile, outfile = tfile, tfile

          # create netcdf version of the input
          Cdo.chainCall("-setname,#{ovar} -copy",:in => @invars[ivar],:out => copyfile,:options => "-f nc")

          # Perform conservative remapping with pregenerated weights
          Cdo.remap([gridfile,weightfile],:in => copyfile,:out => outfile)

          @lock.synchronize { @outvars[ovar] = outfile }
        }
      }
      ths.each {|t| t.join}
    else
      @outvars.each_key {|ovar|
        Dbg.msg("Processing '#{ovar}' ...",@options[:verbose],@options[:debug])

        # determine the grid of the variable
        grid       = @config.selectBy(:outName => ovar).grid[0].to_sym
        gridfile   = @grids[GRID_VARIABLES[grid]]
        weightfile = @weights[GRID_VARIABLES[grid]]

        # remap the variable onto its icon grid provided in the configuration
        ivar              = @config.selectBy(:outName => ovar).inName[0]
        copyfile, outfile = tfile, tfile

        # create netcdf version of the input
        Cdo.chainCall("-setname,#{ovar} -copy",:in => @invars[ivar],:out => copyfile,:options => "-f nc")

        # Perform conservative remapping with pregenerated weights
        Cdo.remap([gridfile,weightfile],:in => copyfile,:out => outfile)

        @outvars[ovar] = outfile
      }
    end
  end

  def verticalInterpolation; end

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

class OmipOptions < PreProcOptions
  def self.parse(args)
    options                        = {}
    options[:verbose]              = false
    options[:openmp]               = 4
    options[:persistent_tempfiles] = false
    options[:threaded]             = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: #{__FILE__} [options]"

      opts.separator ""
      opts.separator "Specific options:"

      # Mandatory arguments
      opts.on("-g", "--grid-file FILE",
              "Use FILE read grid definition/resolution") do |file|
        options[:gridfile] = file
              end
      opts.on("-i", "--input-path PATH",
              "Use PATH for reading input files") do |file|
        options[:inputpath] = file
              end
      opts.on("-o", "--output-file FILE",
              "Use FILE write output variables") do |file|
        options[:outputfile] = file
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
        puts 'omip2icon ' + [0,1].join('.')
        exit
      end

      # development switches: USE AT YOUR OWN RIST!
      opts.on("--persistent-tempfiles",
              "DEVELOPEMENT SWITCH: Create persistent intermediate files instead of cleaning up when finishing the script") {|v|
        options[:persistent_tempfiles] = v
      }

      # Boolean switches
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
      mandatory = [:outputfile, :inputpath, :gridfile]
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

    checkFile(options[:gridfile],'','','')
    unless File.readable?(options[:inputpath])
      warn "Unreadable inputpath '#{options[:inputpath]}!"
      exit
    end

    options
  end
end
class PreProcConfig < ExtCsv

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

  def initialize(mode,datatype,params)
    super(mode,datatype,params)

    # reset column names to internal ones if they differ
    DefaultColumnNames.each {|internalName,externalName|
      next if self.datacolumns.include?(internalName.to_s)
      self.send(internalNama.to_s+'=',self.send(externalName))
      self.delete_field(externalName)
    }
    # check grid IDs
    checkGriddes
    # check and reset the nlevel IDs
    if checkNlevels
      puts "Config Data has only allowed values"
      resetNlevels
    else
      warn "Config Data has non allowed values"
      exit
    end
  end

  def display(verbose=false, debug=false)
    Dbg.msg(self.columns(  :code, :outName, :inName, :info).to_string("tsv",false), verbose, debug)
  end

  def hasCode?(code)
    return true if self.code.include?(code)
    return false
  end
  def hasOutVar?(name)
    return true if self.outName.include?(name)
    return false
  end
  def hasInVar?(name)
    return true if self.inName.include?(name)
    return false
  end

  def gridOfVar(varid)
    case
    when varid.kind_of?(Numeric)
      return self.selectBy(:code => "'#{varid}'").grid[0]
    when varid.kind_of?(String)
      return self.selectBy(:outName => varid).grid[0]
    end
    warn "Could not find grid information for var '#{varid}'!"
    warn "Check config with 'checkConfig'"
    exit
  end

  def hasGridinfo?(varid)
    case
    when varid.kind_of?(Numeric)
      return true if %w[cell edge vertex].include?(conf.selectBy(:code => "'#{varid}'").grid[0])
    when varid.kind_of?(String)
      return true if %w[cell edge vertex].include?(conf.selectBy(:outName => varid).grid[0])
    end
    return false
  end

  def check
    checkGriddes
  end

  def checkGriddes
    return false unless self.grid.uniq.sort == %w[cell edge vertex]
    return true
  end
  def checkNlevels
    check = self.nlevels.uniq - DefaultNlevelIDs.values.flatten
    return false unless check.empty?
    return true
  end
  def resetNlevels
    self.nlevels.each_with_index {|v,i|
      defaultEntry = DefaultNlevelIDs.select {|defk,defv| defv.include?(v)}
      if defaultEntry.empty?
        warn "Could not find valid values for replacing '#{v}'"
        exit
      else
        self.nlevels[i] = defaultEntry.keys[0].to_s
      end
    }
  end
end
class Omip2Icon < PreProc4Icon
  def initialize(options)
    # some preparations to have a single input file
    @lock = Mutex.new

    @options = options
    @options.freeze

    @tempfiles = []

    @invars, @outvars, @grids, @weights = {}, {}, {}, {}

    @ifile, @ipath, @ofile = tfile, @options[:inputpath], @options[:outputfile]

    ifiles = findInputFiles(@ipath)

    # Create single input file
    mergeInputFiles(ifiles)

    # Create internal configuration state
    @config = createConfigFromFile(@ifile,"cell",1)
    pp @config if @options[:debug]

    prepareGridAndWeights

    # Read input variables form File
    @invars = getInputVariables

    # Check which input/output variabels have a valid configuration
    checkConfigOnIOVariables

  end

  def findInputFiles(path)
    ifiles = Dir.glob([path,"*.nc"].join(File::SEPARATOR))
    # skip land see mask files
    ifiles.delete_if {|f| /land_sea/.match(f)}
    ifiles
  end
  def mergeInputFiles(ifiles)
    ofiles = []
    ths = []
    ifiles.each {|ifile|
      ths << Thread.new(ifile) {|_ifile|
        ofile = tfile
        Cdo.settaxis('2001-01-01,12:00:00,1day',:in => _ifile, :out => ofile)
        @lock.synchronize{ ofiles << ofile }
      }
    }
    ths.each {|t| t.join}
    Cdo.merge(:in => ofiles.join(' '), :out => @ifile)
  end
  def createConfigFromFile(ifile,grid,nlevel)
    names = Cdo.showname(:in => ifile)
    configData = {
      :outName => names,
      :inName  => names,
      :grid       => (0...names.size).collect {grid},
      :nlevels     => (0...names.size).collect {nlevel.to_s}
    }
    PreProcConfig.new('hash','plain',configData)
  end
  def verticalInterpolation;end
  def checkConfigOnIOVariables
    @invars.each_key {|k| @outvars[k] = nil}
  end
end

# main script if the file is called directly
if __FILE__ == $0
  options  = OmipOptions.parse(ARGV)
  Cdo.setCdo(options[:cdo]) if options[:cdo]
  pp options if options[:debug]
  Cdo.Debug = options[:debug] if options[:verbose]
  #=======================================================
  p    = Omip2Icon.new(options)
  p.run
end
