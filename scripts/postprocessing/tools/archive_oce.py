#!/usr/bin/env python
import os,sys,glob,shutil,json,subprocess,multiprocessing
from cdo import *
import matplotlib
matplotlib.use('Agg')
import pylab,numpy,dateutil.parser


# ==============================================================================
# USAGE {{{
def usage():
  return """
  # USAGE =====================================================================
  #
  #   ./archive.py EXP=<expname> FILEPATTERN=<pattern> ARCHDIR=<dir> TAG=<rxxxxx>
  #
  # ===========================================================================
  #
  # EXPNAME     : string for naming things
  # TAG         : additional placeholder, e.g. for revisions, if not used, set it to ''
  # FILEPATTERN : quoted wildcard for model result relative to the 'experiments' dir
  # ARCHDIR     : directory, where the archived data is placed (default: ./archive)
  #
  # more options are (keys area llways uppercase)
  # DEBUG       : False/True
  # CALCPSI     : path, where to find the psi calculation and plotting scipt
  # ICONPLOT    : call for icon_plot.ncl, default is 'nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools'
  #
  # #
  # !!!! THIS SCRIPT IS MEANT TO BE STARTED IN THE MODEL RESULTS DIRECTORY !!!!
  """
# }}} ----------------------------------------------------------------------------------
# OPTION HANDLING {{{ =========================================================
def parseOptions():
  options = {
              'GRID'        : 'global',                    # shortcut for the used gridtype(global,box,channel)
              'ARCHDIR'     : './archive',
              'PLOTDIR'     : './plots',
              'DOCTYPE'     : 'pdf',                       # target document format
              'EXP'         : 'oce_mpiom',                 # default experiment name
              'FILEPATTERN' : 'oce_mpiom/oce_mpiom_*.nc*', # default output file pattern
              'DEBUG'       : False,                       # debugging is switched of by default
              'FORCE'       : False,                       # recomputation of verything is switched off by default
                              # the psi processor/plotter
              'CALCPSI'     : '../../scripts/postprocessing/tools/calc_psi.py',
              'TAG'         : 'r1xxxx',                    # addition revision information
#             'ICONPLOT'    : 'nclsh ../../scripts/postprocessing/tools/icon_plot.ncl -altLibDir=../../scripts/postprocessing/tools',
              'ICONPLOT'    : 'nclsh /scratch/mpi/CC/mh0287/users/m300064/builds/remote/icon/gcc/icon-ocean_diagnostics/scripts/postprocessing/tools/icon_plot.ncl -altLibDir=/scratch/mpi/CC/mh0287/users/m300064/builds/remote/icon/gcc/icon-ocean_diagnostics/scripts/postprocessing/tools',
              'PROCS'       : 8,                           # number of threads/procs to be used for parallel operations
              'JOBISRUNNING': True,                        # avoid the last output file/result year by default
              # optional stuff
              'DRYRUN'      : False,                       # with this set to true, the model output is scanned for containing years, only
              'MOCPATTERN'  : 'MOC.*',
#              'MOCPLOTTER'  : '../../scripts/postprocessing/tools/calc_moc.ksh',
              'MOCPLOTTER'  : '/scratch/mpi/CC/mh0287/users/m300064/builds/remote/icon/gcc/icon-ocean_diagnostics/scripts/postprocessing/tools/calc_moc.ksh',
              # options to select special parts od the script
              #'ACTIONS'     : 'archive,preproc,procRegio,plotRegio,plotPsi,plotTf,plotHorz,plotX,plotMoc,plotTSR,finalDoc',
              'ACTIONS'     : 'archive,preproc,plotPsi,plotTf,plotHorz,plotX,plotMoc,plotTSR,finalDoc',
             }

  optsGiven = sys.argv[1:]
  for optVal in optsGiven:
      key,value    = optVal.split('=')
      if key in ['FORCE','DEBUG','JOBISRUNNING','DRYRUN']:
        value = value.lower() in ['true','1']
      if 'PROCS' == key:
        value = int(value)
      if 'ACTIONS' == key:
        value = value.split(',')

      options[key] = value
  return options
# }}} ----------------------------------------------------------------------------------
# HELPER METHODS {{{ =========================================================== 
def dbg(obj):
  if options['DEBUG']:
    print(obj)

""" save internal log """
def dumpLog():
  with open(LOGFILE,"w") as f:
    f.write(json.dumps(LOG))

""" save internal log and exit """
def doExit(value=0):
  # save the restults to omitt reprocessing input files on subsequent calls
  dumpLog()
  sys.exit(value)

""" load the last archiving status ir available """
def loadLog(options):
  if options['FORCE']:
    LOG = {}
    if os.path.exists(LOGFILE):
      os.remove(LOGFILE)
  if os.path.exists(LOGFILE):
    dbg("Load last status from LOGFILE")
    LOG = json.load(open(LOGFILE))
  else:
    dbg("Could not find file logile: '#{@logile}'! Create a new one ....")
    options['FORCE'] = True
    LOG = {}

  return LOG

""" append list of image together """ #{{{
def collectImageToMapByRows(images,columns,ofile):
  imageMap     = {}
  imageColumns = {}
  allTempFiles = []
  for i in range(0,columns):
    imageColumns[i] = []

  for i,image in enumerate(images):
    imageColumns[i%columns].append(image)
  for column in imageColumns:
    _tfile  = tempfile.NamedTemporaryFile(delete=True,prefix='IconArchive')
    allTempFiles.append(_tfile)
    rowFile = _tfile.name+'.png'
    cmd     = 'convert +append %s %s'%(' '.join(imageColumns[column]),rowFile)
    dbg(cmd)
    subprocess.check_call(cmd,shell=True,env=os.environ)
    imageMap[column] = rowFile

  cmd = "convert -append %s %s"%(' '.join(imageMap.values()),ofile)
  dbg(cmd)
  subprocess.check_call(cmd,shell=True,env=os.environ)
  map(lambda x: os.remove(x),imageMap.values()) #}}}

""" warpper around cdo commands for multiprocessing """ #{{{
def showyear(file):
  return cdo.showyear(input = file)[0].split(' ')
def ntime(file):
  return cdo.ntime(input = file)[0]
""" wrapper of 'cdo.yearmean' """
def yearmean(ifile,ofile):
  cdo.yearmean(input = ifile, output = ofile)
# }}}

""" collect the years/ntime, which aree stored in to files """ #{{{
def scanFilesForTheirYears(fileList,procs,log):

  if not log.has_key('yearsOfFiles'):
    log['yearsOfFiles'] = {}

  pool = multiprocessing.Pool(procs)
  results = []


  for ifile in fileList:
    years = pool.apply_async(showyear,[ifile])
    results.append([ifile,years])

  for result in results:
      f = result[0]
      r = result[1]
      log['yearsOfFiles'][f] = r.get()

  pool.close()
  pool.join()

  dumpLog()

def scanFilesForTheirNTime(fileList,procs,log):

  pool = multiprocessing.Pool(procs)
  results = []

  for ifile in fileList:
    ntimes = pool.apply_async(ntime,[ifile])
    results.append([ifile,ntimes])

  for result in results:
      f = result[0]
      r = result[1]
      log[f] = r.get()

  pool.close()
  pool.join()
#}}}
""" create yearly and yearMean filenames """
def getFileNamesForYears(year,archdir,experimentInfo):
  yearFile     = "%s/%s_%s.nc"%(archdir,experimentInfo,year)
  yearMeanFile = "%s/%s_%s_mean.nc"%(archdir,experimentInfo,year)

  return [yearFile, yearMeanFile]

""" input processing for computing regio mean values """
def _computeRegioMean(depth,varname,ifile,vertMask,regioMask,ofile,newData):
  # mask out region
  # vertical interpolation to target depth
  # mean value computaion
  cdo.fldmean(input = '-mul -div -sellevel,%s -selname,%s %s %s %s'%(depth,varname,ifile,vertMask,regioMask),
              output = ofile,
              forceOutput = newData)

""" compute global mean of bias for t s rho """
def globalTempSalRho1D(ifile,initfile,varNames,archdir):
  maskName = 'wet_c'
  ofile    = '/'.join([archdir,'TSR_1D_%s'%(os.path.basename(ifile))])
  cdo.timmean(input = '-fldmean -div -sub -selname,%s %s -selname,%s %s -selname,%s -seltimestep,1 %s'%(varNames,ifile,varNames,initfile,maskName,ifile),
              output = ofile)

""" compute horizontal last 30-year bias for t s rho """
def globalTempSalRho2D(ifile,initfile,varNames,archdir):
  maskName = 'wet_c'
  ofile    = '/'.join([archdir,'TSR_2D_%s'%(os.path.basename(ifile))])
  cdo.timmean(input = '-div -sub -selname,%s -sellevidx,1 %s -selname,%s -sellevidx,1 %s -sellevidx,1 -selname,%s -seltimestep,1 %s'%(varNames,ifile,varNames,initfile,maskName,ifile),
              output = ofile)

""" write a single file for a given years from a list of input files """ #{{{
def grepYear(ifiles,year,archdir,forceOutput,shouldComputeYearMean,experimentInfo):
  yearFile, yearMeanFile = getFileNamesForYears(year,archdir,experimentInfo)

  if (not os.path.exists(yearFile) or forceOutput):
    if (os.path.exists(yearFile)):
      os.remove(yearFile)

    catFiles = []
    for ifile in ifiles:
      ofile = "%s/_catFile_%s_%s"%(archdir,year,os.path.basename(ifile))
      catFiles.append(cdo.selyear(year, input  = ifile, output = ofile))
    if ( 1 == len(ifiles) ):
      cdo.copy(input = catFiles[0], output = yearFile)
    else:
      cdo.cat(input = " ".join(catFiles), output = yearFile)

    if (shouldComputeYearMean):
      cdo.yearmean(input = yearFile,output = yearMeanFile,forceOutput = forceOutput)
    else:
      yearMeanFile = None

    map(lambda x: os.remove(x),catFiles)
  else:
    print("Use existing yearFile '%s'"%(yearFile))
#}}}

""" split input files into yearly files - compute yearmean if desired """
def splitFilesIntoYears(filesOfYears,archdir,forceOutput,expInfo,procs):
  pool      = multiprocessing.Pool(procs)
  yearFiles = []
  for year,files in filesOfYears.iteritems():
    yearFile               = pool.apply_async(grepYear,[files,str(year),archdir,forceOutput,True,expInfo])
    yearFile, yearMeanFile = getFileNamesForYears(str(year),archdir,expInfo)
    yearFiles.append([year,yearFile,yearMeanFile])

  pool.close()
  pool.join()

  return yearFiles

""" compute timmean + fldmean of masked temperature, salinity and potential density """
def computeMaskedTSRMeans1D(ifiles,varList,exp,archdir,procs):
  pool    = multiprocessing.Pool(procs)
  results = []

  initFile  = cdo.selname(','.join(varList),input = '-seltimestep,1 %s'%(ifiles[0]),
                                            output = '/'.join([archdir,'TSR_%s_init'%(exp)]))
  # compute the 1d (vertical) year mean variance to the initial state
  for ifile in ifiles:
    ofile = pool.apply_async(globalTempSalRho1D,[ifile,initFile,','.join(varList),archdir])
    ofile = '/'.join([archdir,'TSR_1D_%s'%(os.path.basename(ifile))])
    results.append(ofile)

  pool.close()
  pool.join()

  merged = '/'.join([archdir,'TSR_1D_%s_complete.nc'%(options['EXP'])])
#TODO if os.path.exists(merged):
#TODO   os.remove(merged)
  merged = cdo.cat(input = ' '.join(sorted(results)),
                   output =  merged)

  if subprocess.check_call("ncrename -d depth_2,depth -v depth_2,depth -O %s"%(merged),shell=True,env=os.environ):
    print("ERROR: ncrename failed")
    print("ERROR: ncrename -d depth_2,depth -v depth_2,depth_2 -O %s"%(merged))
    exit(1)

  return merged

""" compute timmean + fldmean of masked temperature, salinity and potential density """
def computeMaskedTSRMeans2D(ifiles,varList,exp,archdir,procs,log):
# pool    = multiprocessing.Pool(procs)
# results = []

# initFile  = cdo.selname(','.join(varList),input = '-seltimestep,1 %s'%(ifiles[0]),
#                                           output = '/'.join([archdir,'TSR_%s_init'%(exp)]))
#
# # compute the 2d surface temperature + salinity variance to the initial
# # state; relevant are the last 30 years because PHC itselt is a 30-year
# # climatology
# for ifile in ifiles:
#   ofile = pool.apply_async(globalTempSalRho2D,[ifile,initFile,','.join(varList),archdir])
#   ofile = '/'.join([archdir,'TSR_2D_%s'%(os.path.basename(ifile))])
#   results.append(ofile)
#
# pool.close()
# pool.join()
#
# merged = '/'.join([archdir,'TSR_2D_%s_complete.nc'%(options['EXP'])])
# if os.path.exists(merged):
#   os.remove(merged)
# merged = cdo.cat(input = ' '.join(sorted(results)),
#                  output =  merged)
#
# # compute the timmean at the end
# timmeamMerged  = '/'.join([archdir,'TSR_2D_%s_complete_timmean.nc'%(options['EXP'])])
# cdo.timmean(input = merged, output = timmeamMerged)

  # alternative: substract 30year mean values instead of averaging variations
  ofile = '/'.join([archdir,'TSR_2D_%s_complete_timmean.nc'%(options['EXP'])])
  cdo.selname(','.join(varList),
              input  = log['last30YearsMeanBias'],
              output = ofile)
  ofileMasked = '/'.join([archdir,'TSR_2D_%s_complete_timmean_masked.nc'%(options['EXP'])])
  applyMask(ofile,LOG['mask'],ofileMasked)

  return ofileMasked

""" filter for sorting glob results """
def mtime(filename):
  return os.stat(filename).st_mtime

""" create images form 1d txt output """
def plotOnlineDiagnostics(diagnosticFiles,options):
  plots                       = []
  onlineMeanValueFileName     = options['ARCHDIR']+'/TS_last10YearMean.txt'
  onlineMeanValueFile         = open(onlineMeanValueFileName, 'w')
  onlineMeanValueDocumentName = options['PLOTDIR']+'/TS_last10YearMean.pdf'
  # cat all diagnostics together
  diagnosticJoined = options['ARCHDIR']+'/'+'_'.join([options['EXP'],"diagnostics.txt"])
  dbg(diagnosticFiles)

  # read in all values columnwise
  all_diag = []
  with open(diagnosticJoined,'w') as joinedFile:
    for i, diagFile in enumerate(diagnosticFiles):
      f     = open(diagFile,'r')
      input = f.readlines()
      for lineIndex, line in enumerate(input):
        if (0 == lineIndex and 0 != i):
          continue
        newLine          = re.sub('(\s)+', r'\1', line)
        joinedFile.write(newLine)
        newLine          = newLine.rstrip('\r\n')
        input[lineIndex] = newLine
        all_diag.append(newLine.split(' '))
  all_diag = zip(*all_diag)


  header     = open(diagnosticFiles[0],'r').readline().rstrip().split(' ')
  dateIndex  = header.index('date')

  passages = ['gibraltar',
              'denmark_strait',
              'drake_passage',
              'indonesian_throughflow',
              'scotland_iceland',
              'mozambique',
              'framStrait',
              'beringStrait',
              'barentsOpening',
              'agulhas',
              'agulhas_long',
              'agulhas_longer',
             ]
  ice_fields = ['ice_extent_nh','ice_volume_nh']
  columns2Plot = passages + \
                 ice_fields + \
                 ['volume',
                  'total_salinity',
                  'total_temperature',
                  ]
  # filter columns according to what is in the input
  columns2Plot = [ x for x in columns2Plot if x in header ]
  dbg(columns2Plot)
  columns2PlotIndeces = {}
  for column in columns2Plot:
    columns2PlotIndeces[column] = header.index(column)

  dates = all_diag[dateIndex]
  dates = map(dateutil.parser.parse,dates[1:-1])

  # write mean values to txt file
  # start with head row
  leftColumnWith   = max(map(len,columns2Plot))
  rightColumnWidth = 8
  onlineMeanValueFile.write('   '.join(["Passage/Parameter".ljust(leftColumnWith),''.rjust(rightColumnWidth)])+"\n")
  # go on with data  
  for column in columns2Plot:
    data   = all_diag[columns2PlotIndeces[column]]
    values = map(float,numpy.array(data)[1:-1])
    # compute mean value of last 10 years
    last10YearsDates  = []
    last10YearsValues = []
    for i,date in enumerate(dates):
      if str(date.year) in LOG['years'][-11:-1]:
        last10YearsDates.append(date)
        if column in passages:
          last10YearsValues.append(values[i]/1E+9)
        else:
          last10YearsValues.append(values[i])
    if column in passages:
      last10YearMean = str(round(numpy.array(last10YearsValues).mean(),2))
    else:
      last10YearMean = '%.2E' % round(numpy.array(last10YearsValues).mean(),2)

    dbgOutput = ' : '.join([column.ljust(leftColumnWith),
                           last10YearMean.rjust(8)])
    dbg(dbgOutput)
    onlineMeanValueFile.write(dbgOutput+"\n")


    # complete timeseries
    ofile = "%s/TS_%s_%s_%s.png"%(options['PLOTDIR'],column,options['EXP'],options['TAG'])
    if ( not os.path.exists(ofile) or options['FORCE']):
      if column in passages:
        unit = "Sv"
      else:
        unit = ""

      pylab.title("%s , %s : %s TF - all %s years | last 10year mean: %s %s"%(options['EXP'],
                                                                              options['TAG'],
                                                                              column.upper(),
                                                                              len(LOG['years']),
                                                                              str(last10YearMean),
                                                                              unit),
                  fontsize=8)
      pylab.grid()
      if column in ice_fields:
        pylab.plot_date(last10YearsDates, last10YearsValues, linestyle='-',marker=',',linewidth=1.0)  
      else:
        pylab.plot_date(dates, values, linestyle='-',marker=',',linewidth=1.0)  
      pylab.savefig(ofile)
      pylab.clf() 
    plots.append(ofile)
#   # first 10 years
#   ofile = "%s/%s_first10Years_%s_%s.png"%(options['ARCHDIR'],column,options['EXP'],options['TAG'])
#   if ( not os.path.exists(ofile) or options['FORCE']):
#     pylab.title("%s , %s : %s TF - first 10y"%(options['EXP'],options['TAG'],column.upper()),fontsize=8)
#     pylab.grid()
#     pylab.plot_date(dates[0:120], values[0:120], linestyle='-',marker=',',linewidth=1.0)  
#     pylab.savefig(ofile)
#     pylab.clf() 
#   plots.append(ofile)
#   # first 20 years
#   ofile = "%s/%s_first20Years_%s_%s.png"%(options['ARCHDIR'],column,options['EXP'],options['TAG'])
#   if ( not os.path.exists(ofile) or options['FORCE']):
#     pylab.title("%s , %s : %s TF - first 20y"%(options['EXP'],options['TAG'],column.upper()),fontsize=8)
#     pylab.grid()
#     pylab.plot_date(dates[0:240], values[0:240], linestyle='-',marker=',',linewidth=1.0)  
#     pylab.savefig(ofile)
#   plots.append(ofile)

  dbg(plots)
  onlineMeanValueFile.flush()
  onlineMeanValueFile.close()

  # create pdf version of the data table in plots directory 
  # start with postscript
   # cmd = "cat %s |  convert -pointsize 24  label:@-  %s"%(onlineMeanValueFileName,onlineMeanValueDocumentName)
  cmd = "enscript %s -o - | ps2pdf - %s"%(onlineMeanValueFileName,onlineMeanValueDocumentName)
  dbg(cmd)
  if subprocess.check_call(cmd, shell=True, env=os.environ):
    print("Creation of data table %s failed!"%(onlineMeanValueDocumentName))
  if options['DEBUG']:
    print("Created %s"%(onlineMeanValueDocumentName))

  return onlineMeanValueDocumentName

""" process model output for regio mean values """
def processRegionMean(options,regioCodes,regioDepths,regioVars,regioMaskVar):
  regioMeanData = {}
  regioPool     = multiprocessing.Pool(options['PROCS'])
  regioLock     = multiprocessing.Lock()

  # create the regio mask first
  regioMasks    = {}
  for location, regioCode in regioCodes.iteritems():
    ofile     = '/'.join([options['ARCHDIR'],'_'.join(['regioMask',location])+'.nc'])
    ofileTemp = '/'.join([options['ARCHDIR'],'_'.join(['_regioMask',location])+'.nc'])
    cdo.eqc(regioCode,input = '-selname,%s -seltimestep,1 %s'%(regioMaskVar,iFiles[0]),output = ofileTemp)
    regioMasks[location] = cdo.div(input = '%s %s'%(ofileTemp,ofileTemp),output = ofile)
  # create the mask from 3d mask wet_c
  regioVertMasks    = {}
  for depth in regioDepths:
    depth = str(depth)
    ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioVertMask',depth+'m'])+'.nc'])
    regioVertMasks[depth] = cdo.sellevel(depth,input = '-div -selname,wet_c -seltimestep,1 %s -selname,wet_c -seltimestep,1 %s'%(iFiles[0],iFiles[0]),output = ofile)
  # compute the regional mean values
  for location, regioCode in regioCodes.iteritems():
    regioMeanData[location] = {}
    for depth in regioDepths:
      regioMeanData[location][str(depth)] = {}
      for varname in regioVars:
        regioMeanData[location][str(depth)][varname] = {}
        ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,str(depth)+'m'])+'.nc'])
        regioPool.apply_async(_computeRegioMean,[depth,varname,ymFile,regioVertMasks[str(depth)],regioMasks[location],ofile,hasNewFiles])
        regioLock.acquire()
        regioMeanData[location][str(depth)][varname] = ofile 
        regioLock.release()
  regioPool.close()
  regioPool.join()

  return regioMeanData

""" create a single document out of all plots """
def createOutputDocument(plotdir,firstPage,docname,doctype,debug):
  images = glob.glob(plotdir+'/*.png')
  dbg(images)

  ofiles = []
  for image in images:
    ofile = '.'.join([os.path.splitext(image)[0],doctype])
    if subprocess.check_call("convert %s %s"%(image,ofile),shell=True,env=os.environ):
      print("ERROR: Convert of %s to %s failed!"%(image,ofile))
    if debug:
      print("Convert %s %s"%(image,ofile))

    ofiles.append(ofile)

  if ( not '' == firstPage ):
    ofiles.insert(0,firstPage)

  document = '/'.join([plotdir,'.'.join([docname,doctype])])

  if ('ps' == doctype):
    if subprocess.check_call("psmerge -o%s %s"%(document,' '.join(ofiles)),shell=True,env=os.environ):
      print("joining postscripts together into %s failed!")
  elif ('pdf' == doctype):
    if subprocess.check_call("pdftk %s cat output %s"%(' '.join(ofiles),document),shell=True,env=os.environ):
      print("joining pdfs together into %s failed!")
  else:
    print("Unsupported document output type: %s!"%(doctype))
    return

""" consistent yearmean filename """
def yearMeanFileName(archdir,exp):
    return  '/'.join([archdir,'_'.join([exp,'yearmean.nc'])])

""" mask by devision with maskfile """
def applyMask(ifile,maskfile,ofile):
  return cdo.div(input = '%s %s'%(ifile,maskfile), output = ofile)
# }}} --------------------------------------------------------------------------
#=============================================================================== 
# MAIN =========================================================================
# script for archiving (suprise,surpise):
#
# [a] transform output to
#   ymean
#   mommean
#   yealy
# [b] compute diagnostics:
#   psi from 20ym
#   moc from 20ym
#   throughflows from total ym 
#   glob. vertical mixing form total monmean
#
#=============================================================================== 
# INTERNALS {{{ ================================================================
options         = parseOptions()
dbg(options)
LOGFILE         = 'archive.log'
LOG             = loadLog(options)
dbg(LOG)
LOG['options']  = options
cdo             = Cdo()
cdo.debug       = options['DEBUG']
cdo.forceOutput = options['FORCE']
# }}}
# -----------------------------------------------------------------------------
# INPUT HANDLING: {{{
# are the files accessable?
iFiles = glob.glob(options['FILEPATTERN'])
dbg("All available files: "+','.join(iFiles))

# SKIP the last file, if the job is still running
if options['JOBISRUNNING']:
  iFiles.pop()
dbg("Files to be used: "+','.join(iFiles))
if 0 == len(iFiles):
  print usage()
  print("Could not find any result files!")
  doExit(1)
else:
  print("will use \n" + "\n".join(iFiles))

# is the archive/plots dir accessable
for location in ['ARCHDIR','PLOTDIR']:
  if not os.path.isdir(options[location]):
    os.makedirs(options[location])

# look for new files (process only those, unless the FORCE options is not set
if LOG.has_key('yearsOfFiles'):
  processedFiles = LOG['yearsOfFiles'].keys()
else:
  processedFiles = []
newFiles = set(iFiles) ^ set(processedFiles)
hasNewFiles = (len(newFiles) != 0)
if hasNewFiles:
  dbg("New Files found:")
  dbg(newFiles)

# get the data dir from iFiles for later use
LOG['dataDir'] = os.path.abspath(os.path.dirname(iFiles[0]))
# }}}
# =======================================================================================
# DATA SPLITTING {{{ ====================================================================
if options['FORCE']:
  scanFilesForTheirYears(iFiles,options['PROCS'],LOG)
else:
  scanFilesForTheirYears(newFiles,options['PROCS'],LOG)
dbg(LOG['yearsOfFiles'])

# compute all contributing simulation years
allYears = set(); _allYears = LOG['yearsOfFiles'].values()
for years in _allYears:
  allYears.update(years)
allYears = list(allYears)
allYears.sort()
# drop the last one if job is running
if options['JOBISRUNNING']:
  allYears.pop()
dbg(allYears)
# collect all files, which contain data for given years
filesOfYears = {}
for year in allYears:
  yearFiles = []
  for _file, _years in LOG['yearsOfFiles'].iteritems():
    if year in _years:
      yearFiles.append(_file)
  dbg(yearFiles)
  filesOfYears[year] = yearFiles
LOG['filesOfYears'] = filesOfYears

if options['DRYRUN']:
  doExit()

# process each year separately: collect yearly data, compute year mean files
splitInfo = splitFilesIntoYears(LOG['filesOfYears'],
                                options['ARCHDIR'],
                                options['FORCE'],
                                options['EXP'],
                                options['PROCS'])
years, yearMeanFiles = [],[]
LOG['meansOfYears'] = {}
for _info in splitInfo:
  year,yearlyFile,yearMeanFile = _info[0], _info[1], _info[2]
  LOG[year] = yearlyFile
  years.append(year)
  yearMeanFiles.append(yearMeanFile)
  LOG['meansOfYears'][year] = yearMeanFile

years.sort()
yearMeanFiles.sort()
LOG['years']      = years
LOG['splityear?'] = True
dumpLog()
dbg(LOG)
# }}} ===================================================================================
# COMPUTE INITIAL VALUES FILES FOR LATER BIASES {{{ =====================================
LOG['init'] = cdo.seltimestep(1,input =  iFiles[0], output = '/'.join([options['ARCHDIR'],'%s_init'%(options['EXP'])]))
# }}} ===================================================================================
# COMPUTE CELL MASK FOR LATER APPLICATION {{{ ===========================================
LOG['mask'] = cdo.selname('wet_c',input = '-seltimestep,1 %s'%(iFiles[0]), output = '/'.join([options['ARCHDIR'],'%s_mask.nc'%(options['EXP'])]))
# }}} ===================================================================================
# COMPUTE SINGLE YEARMEAN FILES {{{ =====================================================
ymFile = yearMeanFileName(options['ARCHDIR'],options['EXP'])
if ( not os.path.exists(ymFile) or options['FORCE'] or hasNewFiles):
  if (os.path.exists(ymFile)):
    os.remove(ymFile)
  cdo.cat(input=" ".join(yearMeanFiles),output=ymFile)
  # rm ymFiles
  #map(lambda x: os.remove(x),ymFiles)
else:
  print("Use existing ymFile: "+ymFile)
# }}} ===================================================================================
# COMPUTE SINGLE MEAN FILE FROM THE LAST COMPLETE 30 YEARS {{{ =====================================================
LOG['last30YearsMean']     = cdo.timmean(input = '-selyear,%s %s'%(','.join(LOG['years'][-32:-2]),yearMeanFileName(options['ARCHDIR'],options['EXP'])),
                                     output = '%s/last30YearsMean_%s_%s-%s.nc'%(options['ARCHDIR'],options['EXP'],LOG['years'][-32],LOG['years'][-2]))
LOG['last30YearsMeanBias'] = cdo.sub(input = '-timmean -selyear,%s %s -seltimestep,1 %s'%(','.join(LOG['years'][-32:-2]),yearMeanFileName(options['ARCHDIR'],options['EXP']),LOG['init']),
                                     output = '%s/last30YearsMeanBias_%s_%s-%s.nc'%(options['ARCHDIR'],options['EXP'],LOG['years'][-32],LOG['years'][-2]))
# }}} ===================================================================================
# PREPARE INPUT FOR PSI CALC {{{
# collect the last 20 years if there are more than 40 years, last 10 otherwise
if len(LOG['years']) > 40:
  nyears4psi = 20
else:
  nyears4psi = 10
years4Psi = LOG['years'][-(nyears4psi+2):-1]
dbg(LOG['years'])
dbg(years4Psi)
yearInfo  = '-'.join([years4Psi[0],years4Psi[-1]])
uvintName = "u_vint_acc"
uvintFile = '/'.join([options['ARCHDIR'],'_'.join([uvintName,yearInfo])+'.nc'])
cdo.timmean(input = "-selname,%s -selyear,%s/%s %s"%(uvintName,years4Psi[0],years4Psi[-1],ymFile),
            output = uvintFile)
# }}} ----------------------------------------------------------------------------------
# PREPARE INPUT FOR PROFILES {{{
# TODO: SWITCHED OFF
#varNames = ['t_acc','s_acc','u_acc','v_acc']
##  collect the data first (first 10 years)
#varFiles = []
#for year in LOG['years'][0:10]:
## varFiles.append(cdo.selname(' '.join(varNames),
##   input=" -remapnn,lon=-30.0_lat=-65.0 %s"%(LOG[year]),
##   output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))
#  varFiles.append(cdo.remapnn('lon=-30.0_lat=-65.0' , 
#                              options = '-O',
#                              input = ' -selname,%s %s'%(','.join(varNames),LOG[year]),
#                              output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))
#
#varData = cdo.cat(input=' '.join(varFiles),
#                    output="%s/_data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP']))
## add absolute velocity
#velocityData = cdo.expr("'velocity=sqrt(u_acc*u_acc+v_acc*v_acc);'",input=varData,output='/'.join([options['ARCHDIR'],
#                                                                                                   'vel4profile.nc']))
#profileOutputFile = "%s/data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP'])
#if os.path.exists(profileOutputFile):
#  os.remove(profileOutputFile)
#profileData = cdo.merge(input='%s %s'%(varData,velocityData),output=profileOutputFile)
#if subprocess.check_call("ncrename -d depth_2,depth -O %s"%(profileOutputFile),shell=True,env=os.environ):
#  print("ERROR: ncrename failed")
# }}} ----------------------------------------------------------------------------------
# PREPARE INPUT FOR REGIONAL MEANS from yearMean output {{{
# for global grid only
if 'procRegio' in options['ACTIONS']:
  # setup
  regioCodes    = {'NorthAtlantic' : 4,'TropicalAtlantic' : 5, 'SouthernOcean' : 6}
  regioDepths   = [110,215,895,2200]
  regioVars     = ['t_acc','s_acc']
  regioMaskVar  = 'regio_c'
  regioMeanData = processRegionMean(options,regioCodes,regioDepths,regioVars,regioMaskVar)
# }}} ----------------------------------------------------------------------------------
# PREPARE INPUT FOR MOC PLOT {{{
# collect all MOC files
dbg(options['MOCPATTERN'])
mocFiles        = sorted(glob.glob(options['MOCPATTERN']),key = mtime)
# default is to take the mean value ofthe at 10 years as input for the plotscript
# this means 120 months, with monthly output, this is 120 timesteps
mocNeededNSteps = 120
mocLog          = {}
scanFilesForTheirNTime(mocFiles,options['PROCS'],mocLog)
dbg(mocLog)
# check for the numbe rof timesteps in the last moc file
mocLastNtime    = int(mocLog[mocFiles[-1]])
mocMeanFile     = '/'.join([options['ARCHDIR'],'mocMean'])
if mocNeededNSteps <= mocLastNtime:
  # take the last 120 values for timmeaninput
  mocMeanFile = cdo.timmean(input = "-seltimestep,%s/%s %s"%(mocLastNtime-mocNeededNSteps+1,mocLastNtime,mocFiles[-1]),
                            output = mocMeanFile)
else:
  mocMeanFile = cdo.timmean(input = mocFiles[-1], output = mocMeanFile)
dbg(mocMeanFile)
# }}} -----------------------------------------------------------------------------------
# PREPARE DATA FOR T,S,RHOPOT BIAS TO INIT {{{ ------------------------------------------
# target is a year mean files of fldmean data, but mean value computation
# should come at the very end of the processing chain
t_s_rho_Input_1D = []
t_s_rho_Input_2D = []
for year in LOG['years']:
  t_s_rho_Input_1D.append(LOG[year])
t_s_rho_Output_1D = computeMaskedTSRMeans1D(t_s_rho_Input_1D,['t_acc','s_acc','rhopot_acc'],
                                            options['EXP'],options['ARCHDIR'],options['PROCS'])
dbg(t_s_rho_Output_1D)
for year in LOG['years'][-32:-2]:
  t_s_rho_Input_2D.append(LOG[year])
t_s_rho_Output_2D = computeMaskedTSRMeans2D(t_s_rho_Input_2D,['t_acc','s_acc','rhopot_acc'],
                                            options['EXP'],options['ARCHDIR'],options['PROCS'],LOG)
dbg(t_s_rho_Output_2D)
# }}} -----------------------------------------------------------------------------------
# DIAGNOSTICS ===========================================================================
# PSI {{{
plotFile = options['PLOTDIR']+'/'+"_".join(["psi",yearInfo,options['EXP'],options['TAG']+'.png'])
if not os.path.exists(plotFile):
  cmd = '%s %s %s'%(options['CALCPSI'], uvintFile, "LEVELS=20 AREA=%s PLOT=%s"%(options['GRID'],plotFile))
  dbg(cmd)
  if subprocess.check_call(cmd,shell=True,env=os.environ):
    print("ERROR: CALCPSI failed")
# }}} ----------------------------------------------------------------------------------
# HORIZONTAL PLOTS: t,s,u,v,abs(velocity) {{{
horizontalConfig = {
  'varNames'      : ['t_acc','s_acc','h_acc','u_acc','v_acc'],
  'iFile'         : LOG[LOG['years'][-2]],              # use the last COMPLETE year #TODO: dont show monthly mean!!!
  'availableVars' : cdo.showname(input = LOG[LOG['years'][-2]])[0].split(' '),
  'sizeOpt'       : '-xsize=1200 -ysize=800',
}
for varname in horizontalConfig['varNames']:
  if ( varname in horizontalConfig['availableVars'] ):
    # surface plot
    oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-atSurface_'+'_'.join([options['EXP'],options['TAG']])])
    if ( not os.path.exists(oFile+'.png') or options['FORCE'] or hasNewFiles ):
      title = '%s: last yearmean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(horizontalConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType=png',horizontalConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
      elif ('chanel' == options['GRID']):
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')
      else:
        cmd.append('')
      dbg(' '.join(cmd))
      subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)

    # plot for roughly 100m depth
    #   skip sea srface height
    if ('h_acc' == varname ):
      continue

    oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-at100m_'+'_'.join([options['EXP'],options['TAG']])])
    if ( not os.path.exists(oFile+'.png') or options['FORCE'] or hasNewFiles ):
      title = '%s: last yearmean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(horizontalConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType=png',horizontalConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
        cmd.append('-levIndex=5')
      elif ('global' == options['GRID']):
        cmd.append('-levIndex=4')
      else:
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')
      dbg(' '.join(cmd))
      subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)

# }}} ----------------------------------------------------------------------------------
# THROUGH FLOWS {{{
# for global grid only
if ( 'global' == options['GRID'] ):
  diagnosticFiles = sorted(glob.glob(os.path.sep.join([LOG['dataDir'],"oce_diagnostics-*txt"])),key=mtime)
  diagnosticTable = plotOnlineDiagnostics(diagnosticFiles,options)
else:
  diagnosticTable = ''
# }}} ----------------------------------------------------------------------------------
# ATLANTIC X-Section: t,s,rhopot  {{{ ================================
# for global grid only
if ( 'global' == options['GRID'] ):
  for varname in ['t_acc','s_acc','rhopot_acc']:
    # plot last yearmean state
    oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_'+'_'.join([options['EXP'],options['TAG']])])
    iFile4XSection = LOG['meansOfYears'][LOG['years'][-2]] # take last yearmean file from a complete year
    # mask by devision
    iFile4XSection = cdo.div(input  = '-selname,%s %s %s'%(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                             output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
    if ( 'rhopot_acc' == varname ):
      # substract 1000
      iFile4XSection  = cdo.subc(1000.0,
                                 input = '-selname,%s %s'%(varname,iFile4XSection),
                                 output = os.path.splitext(iFile4XSection)[0]+'_%s_subc1000.nc'%(varname))
    if ( not os.path.exists(oFile+'.png') or options['FORCE']):
      title = '%s: last yearmean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType=png',
             '-resolution=r180x90',
             '-selPoints=150',
             '-rStrg="-"',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      dbg(' '.join(cmd))
      subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)

    # plot bias to initialization
    oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_BiasToInit'+'_'.join([options['EXP'],options['TAG']])])
    iFile4XSection = LOG['last30YearsMeanBias'] # take last 30 yearmean bias to init
    # mask by devision
    iFile4XSection = cdo.div(input  = '-selname,%s %s %s'%(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                             output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
    if ( not os.path.exists(oFile+'.png') or options['FORCE']):
      title = '%s: last yearmean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType=png',
             '-resolution=r180x90',
             '-selPoints=150',
             '-rStrg="-"',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      dbg(' '.join(cmd))
      subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)
# }}} ----------------------------------------------------------------------------------
# SOUTH OCEAN t,s,y,v profile at 30w, 65s  {{{ ================================
#  create hovmoeller-like plots
#for varname in ['t_acc','s_acc','u_acc','v_acc','velocity']:
#  # run icon_plot.ncl
#  oFile = '/'.join([options['ARCHDIR'],varname+'_profileAt-30W-65S'+'_'.join([options['EXP'],options['TAG']])])
#  iFile4XSection = profileData
#  if ( not os.path.exists(oFile+'.png') or options['FORCE']):
#    title = '%s: 10y prfile at 30W,65S '%(options['EXP'])
#    cmd = [options['ICONPLOT'],
#           '-iFile=%s'%(iFile4XSection),
#           '-hov',
#           '-varName=%s'%(varname),
#           '-oType=png',
#           '-rStrg="-"',
#           '-withLineLabels',
#           '-tStrg="%s"'%(title),
#           '-oFile=%s'%(oFile)]
#    dbg(' '.join(cmd))
#    subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)
# }}} ----------------------------------------------------------------------------------
# REGIO MEAN PROFILES {{{ ================================
# for global grids only
if ( 'procRegio' in options['ACTIONS'] and 'plotRegio' in options['ACTIONS']):
  regioPlotNames = {'t_acc' : 'Temperature','s_acc' : 'Salinity'}
  for location, regioCode in regioCodes.iteritems():
    for varname in regioVars:
      imageCollection = []
      for depth in regioDepths:
        ifile = regioMeanData[location][str(depth)][varname]

        data  = cdo.readArray(ifile,varname).flatten()

        dates = re.sub('(\s)+', r'\1', cdo.showdate(input = ifile)[0]).split(' ')
        dates = map(dateutil.parser.parse,dates)

        title = '%s at %s: %s [%s,%s]'%(location,str(depth)+'m',varname,options['EXP'],options['TAG'])
        ofile = '/'.join([options['ARCHDIR'],'_'.join(['.regioMean',location,varname,str(depth)+'m',options['EXP'],options['TAG']])+'.png'])

        unit  = cdo.showunit(input = ifile)[0]

        pylab.title(title,fontsize=9)
        pylab.grid()
        pylab.xlabel("Years")
        pylab.ylabel('%s [%s]'%(regioPlotNames[varname],unit))
        pylab.plot_date(dates, data, linestyle='-',marker='.')  
        pylab.savefig(ofile)
        imageCollection.append(ofile)
        pylab.clf()
      collectImageToMapByRows(imageCollection,
                              2,
                              '/'.join([options['PLOTDIR'],'_'.join(['regioMean',location,varname,options['EXP'],options['TAG']+'.png'])])) 
# }}} ----------------------------------------------------------------------------------
# MOC PLOT {{{
mocPlotCmd = '%s %s'%(options['MOCPLOTTER'],mocMeanFile)
mocPlotSetup = {
  'TITLE' : 'ICON: %s : MOC : %s'%(options['EXP'],mocFiles[-1]),
  'IFILE' : mocMeanFile,
  'OFILE' : '/'.join([options['PLOTDIR'],'MOC_%s'%(options['EXP'])]),
  'OTYPE' : 'png',
}
for k,v in mocPlotSetup.iteritems():
  os.environ[k] = v
if subprocess.check_call(mocPlotCmd,shell=True,env=os.environ):
  print("ERROR: MOCPLOT failed")
# environment cleanup
for k in mocPlotSetup.keys():
  os.environ.pop(k)
# }}} ----------------------------------------------------------------------------------
# T S RHOPOT BIAS PLOT {{{
# global mean bias over depth and time
t_s_rho_PlotSetup = {
  't_acc'      : {'maxVar' : '3.0', 'minVar' : '-3.0' , 'numLevs' : '20'},
  's_acc'      : {'maxVar' : '0.2', 'minVar' : '-0.2' , 'numLevs' : '16'},
  'rhopot_acc' : {'maxVar' : '0.6', 'minVar' : '-0.6' , 'numLevs' : '24'},
}
for varname in t_s_rho_PlotSetup.keys():
  oFile = '/'.join([options['PLOTDIR'],varname+'_biasToInit_inDepth_overTime'+'_'.join([options['EXP'],options['TAG']])])
  #if ( not os.path.exists(oFile+'.png') or options['FORCE'] or hasNewFiles ):
  title = '%s: %s bias to init '%(options['EXP'],varname)
  cmd = [options['ICONPLOT'],
         '-iFile=%s'%(t_s_rho_Output_1D),
         '-varName=%s'%(varname),
         '-oType=png',
         '-hov',
         '-rStrg="-"',
         '-bStrg="-"',
         '-lStrg=%s'%(varname),
#        '-withLineLabels',
         '+withLines',
         '-minVar=%s'%(t_s_rho_PlotSetup[varname]['minVar']),
         '-maxVar=%s'%(t_s_rho_PlotSetup[varname]['maxVar']),
         '-numLevs=%s'%(t_s_rho_PlotSetup[varname]['numLevs']),
         '-colormap=BlWhRe',
  #      '-tStrg="%s"'%(title),
         '-oFile=%s'%(oFile)]
  print(cmd)
  dbg(' '.join(cmd))
  if subprocess.check_call(' '.join(cmd),shell=True,env=os.environ):
    print('CMD: %s has failed!')

# surface bias of the last 30 year mean
horizontalConfig = {
  'varNames'      : ['t_acc','s_acc'],
  'iFile'         : t_s_rho_Output_2D,
  'availableVars' : cdo.showname(input = t_s_rho_Output_2D)[0].split(' '),
  'sizeOpt'       : '-xsize=1200 -ysize=800',
}
for varname in horizontalConfig['varNames']:
  if ( varname in horizontalConfig['availableVars'] ):
    # surface plot
    oFile = '/'.join([options['PLOTDIR'],varname+'_biasToInit_atSurface_last30YearsMean'+'_'.join([options['EXP'],options['TAG']])])
    if ( not os.path.exists(oFile+'.png') or options['FORCE'] or hasNewFiles ):
      title = '%s: 30 year mean surface bias to init'%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-isIcon',
             '-iFile=%s'%(horizontalConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType=png',horizontalConfig['sizeOpt'],
             '-rStrg="-"',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
      elif ('chanel' == options['GRID']):
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')
      else:
        cmd.append('')
      dbg(' '.join(cmd))
      subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)
# }}} ----------------------------------------------------------------------------------
# FINAL DOCUMENT CREATION {{{ ===========================================================
createOutputDocument(options['PLOTDIR'],diagnosticTable,'_'.join(['ALL',options['EXP'],options['TAG']]),options['DOCTYPE'],options['DEBUG'])
# }}} ----------------------------------------------------------------------------------
#
#
# vim:fdm=marker
