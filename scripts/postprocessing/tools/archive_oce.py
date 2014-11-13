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
              'OFORMAT'     : 'png',                       # target document format
              'EXP'         : 'oce_mpiom',                 # default experiment name
              'FILEPATTERN' : 'oce_mpiom/oce_mpiom_*.nc*', # default output file pattern
              'DEBUG'       : False,                       # debugging is switched of by default
              'FORCE'       : False,                       # recomputation of verything is switched off by default
                              # the psi processor/plotter
              'CALCPSI'     : '../../scripts/postprocessing/tools/calc_psi.py',
              'TAG'         : 'r1xxxx',                    # addition revision information
              'ICONPLOT'    : 'nclsh ../../scripts/postprocessing/tools/icon_plot.ncl -altLibDir=../../scripts/postprocessing/tools -remapOperator=remapycon',
              'PROCS'       : 8,                           # number of threads/procs to be used for parallel operations
              'JOBISRUNNING': True,                        # avoid the last output file/result year by default
              # optional stuff
              'DRYRUN'      : False,                       # with this set to true, the model output is scanned for containing years, only
              'MOCPATTERN'  : 'MOC.*',
              'MOCPLOTTER'  : '../../scripts/postprocessing/tools/calc_moc.ksh',
              # options to select special parts od the script
              'ACTIONS'     : 'archive,preproc,procMoc,plotMoc,procRegio,plotRegio,plotTf,plotHorz,plotX,procTSR,plotTSR,plotPsi,procIce,plotIce',#finalDoc',
              #'ACTIONS'     : 'archive,preproc,procRegio,plotRegio,plotTf,plotHorz,plotX,plotTSR',#plotPsi',#finalDoc',
              #'ACTIONS'     : 'archive,preproc,plotHorz',#finalDoc',
#             'ACTIONS'     : 'archive,preproc,plotPsi,plotTf,plotHorz,plotX,plotMoc,plotTSR,finalDoc',
             }

  optsGiven = sys.argv[1:]
  for optVal in optsGiven:
    if ( 0 <= optVal.find('=') ):
      if (optVal.find('=') == optVal.rfind('=')):
        key,value = optVal.split('=')
      else:
        firstSeperatorIndex = optVal.find('=')
        key   = optVal[:(optVal.find('='))]
        value = optVal[firstSeperatorIndex+1:]

      if key in ['FORCE','DEBUG','JOBISRUNNING','DRYRUN']:
        value = value.lower() in ['true','1']
      if 'PROCS' == key:
        value = int(value)
      if 'ACTIONS' == key:
        value = value.split(',')

      options[key] = value
    else:
      print("Cannot parse options '%s'!"%(optVal))
      sys.exit(1)

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
def _computeRegioMean(depth,varname,ifile,myMask,ofile,newData):
  # mask out region
  # vertical interpolation to target depth
  # mean value computaion
  cdo.fldmean(input = '-div -sellevel,%s -selname,%s %s %s'%(depth,varname,ifile,myMask),
              output = ofile,
              forceOutput = newData)

""" compute global mean of bias for t s rho """
def globalTempSalRho1D(ifile,initfile,maskfile,varNames,ofile):
  cdo.timmean(input = '-fldmean -div -sub -selname,%s %s %s %s'%(varNames,ifile,initfile,maskfile),
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
def computeMaskedTSRMeans1D(ifiles,varList,initFile,maskFile,exp,archdir,procs):
  pool    = multiprocessing.Pool(procs)
  results = []

  initFile  = cdo.selname(','.join(varList),input = initFile, output = '/'.join([archdir,'TSR_%s_init'%(exp)]))
  # compute the 1d (vertical) year mean variance to the initial state
  for ifile in ifiles:
    _ofile = '/'.join([archdir,'TSR_1D_%s'%(os.path.basename(ifile))])
    ofile  = pool.apply_async(globalTempSalRho1D,[ifile,initFile,maskFile,','.join(varList),_ofile])
    ofile  = _ofile
    results.append(ofile)

  pool.close()
  pool.join()

  merged = '/'.join([archdir,'TSR_1D_%s_complete.nc'%(options['EXP'])])
  if os.path.exists(merged):
    os.remove(merged)
  merged = cdo.cat(input = ' '.join(sorted(results)),
                   output =  merged)

  # rename vertical axis if present - ignore exit status (return code differs depending on the ncp version)
  if subprocess.call("ncrename -d depth_2,depth -v depth_2,depth -O %s"%(merged),shell=True,env=os.environ):
    print("ERROR: ncrename failed")
    print("ERROR: ncrename -d depth_2,depth -v depth_2,depth_2 -O %s"%(merged))
#   exit(1)

  return merged

""" compute timmean + fldmean of masked temperature, salinity and potential density """
def computeMaskedTSRMeans2D(ifiles,varList,exp,archdir,procs,log):
  # alternative: substract 30year mean values instead of averaging variations
  ofile = '/'.join([archdir,'TSR_2D_%s_complete_timmean.nc'%(options['EXP'])])
  cdo.selname(','.join(varList),
              input  = log['last30YearsMeanBias'],
              output = ofile)
  ofileMasked = '/'.join([archdir,'TSR_2D_%s_complete_timmean_masked.nc'%(options['EXP'])])
  applyMask(ofile,LOG['mask'],ofileMasked)

  return ofileMasked
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


""" filter for sorting glob results """
def mtime(filename):
  return os.stat(filename).st_mtime

""" split filename into dirname, basenameWithoutExtension, extension """
def splitFileName(filename):
  dirname                     = os.path.dirname(filename)
  filenameWithoutExt, extname = os.path.splitext(os.path.basename(filename))

  return [dirname,filenameWithoutExt,extname]

""" wrapper for single command line execution """
def call(command):
  proc = subprocess.Popen(command,
      shell  = True,
      stderr = subprocess.PIPE,
      stdout = subprocess.PIPE)
  retvals = proc.communicate()
  print(retvals)

""" execute list of system command in parallel """
def executeInParallel(commandList,procs):

  if not commandList:
    return

  pool = multiprocessing.Pool(procs)
  results = []

  for command in commandList:
    pool.apply_async(call,[command])

  pool.close()
  pool.join()

""" create images form 1d txt output """
def plotOnlineDiagnostics(diagnosticFiles,options):
  plots                       = []
  onlineMeanValueFileName     = options['ARCHDIR']+'/TS_last10YearMean_%s.txt'%(options['EXP'])
  onlineMeanValueFile         = open(onlineMeanValueFileName, 'w')
  onlineMeanValueDocumentName = options['PLOTDIR']+'/TS_last10YearMean_%s.pdf'%(options['EXP'])
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
  lastCompleteYear = dates[-1].year - 1
  dbg("Online diagnostics:")
  dbg("lastCompleteYear:")
  dbg(lastCompleteYear)

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
      if (lastCompleteYear - 10 < date.year) and (date.year <= lastCompleteYear):
        #dbg(date)
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
    #if ( not os.path.exists(ofile) or options['FORCE']):
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
def processRegionMean(options,mask3D,init,yearMeanFile,regioCodes,regioDepths,regioVars,regioMaskVar):
  regioMeanData = {}
  regioPool     = multiprocessing.Pool(options['PROCS'])
  regioLock     = multiprocessing.Lock()

  # create the regio mask first
  regioMasks    = {}
  for location, regioCode in regioCodes.iteritems():
    ofile     = '/'.join([options['ARCHDIR'],'_'.join(['regioMask',location])+'.nc'])
    ofileTemp = '/'.join([options['ARCHDIR'],'_'.join(['_regioMask',location])+'.nc'])
    cdo.eqc(regioCode,input = '-selname,%s %s'%(regioMaskVar,init),output = ofileTemp)
    regioMasks[location] = ofileTemp#cdo.div(input = '%s %s'%(ofileTemp,ofileTemp),output = ofile)

  # create the mask from 3d mask wet_c
  regioVertMasks    = {}
  for depth in regioDepths:
    depth = str(depth)
    ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioVertMask',depth+'m'])+'.nc'])
    regioVertMasks[depth] = cdo.sellevel(depth,input = mask3D,output = ofile)

  # compute the regional mean values
  for location, regioCode in regioCodes.iteritems():
    regioMeanData[location] = {}
    for depth in regioDepths:
      regioMeanData[location][str(depth)] = {}

      # create a single mask file out of regional and vertical mask
      myMask = cdo.mul(input=' '.join([regioMasks[location],regioVertMasks[str(depth)]]),
                       output='/'.join([options['ARCHDIR'],'_'.join(['regioMask',location,str(depth)+'m'])+'.nc']))

      for varname in regioVars:
        regioMeanData[location][str(depth)][varname] = {}
        ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,str(depth)+'m'])+'.nc'])
        regioPool.apply_async(_computeRegioMean,[depth,varname,yearMeanFile,myMask,ofile,hasNewFiles])
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
def yearMeanFileName(archdir,exp,firstYear,lastYear):
    return  '/'.join([archdir,'_'.join([exp,'%s-%s'%(firstYear,lastYear),'yearmean.nc'])])

def yearMonMeanFileName(archdir,exp,firstYear,lastYear):
    return  '/'.join([archdir,'_'.join([exp,'%s-%s'%(firstYear,lastYear),'yearmonmean.nc'])])

""" mask by devision with maskfile """
def applyMask(ifile,maskfile,ofile):
  return cdo.div(input = '%s %s'%(ifile,maskfile), output = ofile)

""" X-Section plot """
def plotXSection(plotConfig,options):
  cmd = [options['ICONPLOT']]
  for key,value in plotConfig.iteritems():
      cmd.append('%s=%s'%(key,value))
#       '-iFile=%s'%(t_s_rho_Output_1D),
#       '-varName=%s'%(varname),
#       '-oType='+options['OFORMAT'],
#       '-hov',
#       '-rStrg="-"',
#       '-bStrg="-"',
#       '-lStrg=%s'%(varname),
#       '-withLineLabels',
#       '+withLines',
#       '-minVar=%s'%(t_s_rho_PlotSetup[varname]['minVar']),
#       '-maxVar=%s'%(t_s_rho_PlotSetup[varname]['maxVar']),
#       '-numLevs=%s'%(t_s_rho_PlotSetup[varname]['numLevs']),
#       '-colormap=BlWhRe',
##      '-tStrg="%s"'%(title),
#       '-oFile=%s'%(oFile)]
#print(cmd)
#dbg(' '.join(cmd))
#if subprocess.check_call(' '.join(cmd),shell=True,env=os.environ):
#  print('CMD: %s has failed!')

def plotHorizontal(plotConfig,options,hasNewFiles):
  plots = []
  for varname in plotConfig['varNames']:
    _dirname, _filename, _extname = splitFileName(plotConfig['iFile'])
    if ( varname in plotConfig['availableVars'] ):
      # surface plot
      oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-atSurface_'+'_'.join([_filename,plotConfig['tag']])])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(plotConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],plotConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="%s"'%(plotConfig['title']),
             '-oFile=%s'%(oFile)]

      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
      elif ('chanel' == options['GRID']):
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')
      else:
        cmd.append('')

      if ('maskFile' in plotConfig):
        cmd.append('-maskFile=%s'%(plotConfig['maskFile']))

      if ('mapType' in plotConfig):
        cmd.append('-mapType=%s'%(plotConfig['mapType']))

      if ('colormap' in plotConfig):
        cmd.append('-colormap=%s'%(plotConfig['colormap']))

      if ('opts' in plotConfig):
        cmd.append(plotConfig['opts'])

      if ('limits' in plotConfig):
        limits = plotConfig['limits'][varname]
        if 'plotLevs' in limits:
          cmd.append('-plotLevs=%s '%(limits['plotLevs']))
        else:
          cmd.append('-minVar=%s -maxVar=%s %s'%(limits['minVar'],limits['maxVar'],limits['mode']))

      cmd = ' '.join(cmd)
      dbg(cmd)
      plots.append(cmd)

      # plot for roughly 30m,100m and 200m depth
      #   skip sea surface height
      if ( 1 == len(cdo.showlevel(input = "-selname,%s -seltimestep,1 %s"%(varname,plotConfig['iFile'])))):
        continue

      oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-at100m_'+'_'.join([_filename,plotConfig['tag']])])
      title = '%s: last year mean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(plotConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],plotConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="%s"'%(plotConfig['title']),
             '-oFile=%s'%(oFile)]

      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
        cmd.append('-levIndex=5')
      elif ('global' == options['GRID']):
        cmd.append('-levIndex=8')
      else:
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')

      if ('maskFile' in plotConfig):
        cmd.append('-maskFile=%s'%(plotConfig['maskFile']))

      if ('mapType' in plotConfig):
        cmd.append('-mapType=%s'%(plotConfig['mapType']))

      if ('colormap' in plotConfig):
        cmd.append('-colormap=%s'%(plotConfig['colormap']))

      if ('opts' in plotConfig):
        cmd.append(plotConfig['opts'])

      if ('limits' in plotConfig):
        limits = plotConfig['limits'][varname]
        if 'plotLevs' in limits:
          cmd.append('-plotLevs=%s '%(limits['plotLevs']))
        else:
          cmd.append('-minVar=%s -maxVar=%s %s'%(limits['minVar'],limits['maxVar'],limits['mode']))

      cmd = ' '.join(cmd)
      dbg(cmd)
      plots.append(cmd)

  executeInParallel(plots,options['PROCS'])
  return
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
plotCommands    = []
cdo             = Cdo()
cdo.cdfMod = 'netcdf4'
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

#------------------------------------------------------------------------------
# BASIC PLOTSETUP FOR MAIN VARIABLES
PlotConfig =  {
  't_acc'      : {'plotLevs' : '-2,-1,-0,1,2,5,10,15,20,25,30'},
  's_acc'      : {'plotLevs' : '20,25,28,30,32,33,34,34.5,35,35.5,36,36.5,37,38,40'},
  'rhopot_acc' : {'plotLevs' : '20,25,28,30,32,34,36,38,40'},
  'u_acc'      : {'plotLevs' : '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  'v_acc'      : {'plotLevs' : '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  'h_acc'      : {'plotLevs': '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  }
PlotConfigBias =  {
  't_acc'      : {'maxVar' : '3.0', 'minVar' : '-3.0' , 'numLevs' : '20'},
  's_acc'      : {'maxVar' : '0.2', 'minVar' : '-0.2' , 'numLevs' : '16'},
  'rhopot_acc' : {'maxVar' : '0.6', 'minVar' : '-0.6' , 'numLevs' : '24'},
  }
# }}}
# =======================================================================================
# DATA SPLITTING {{{ ====================================================================
if 'preproc' in options['ACTIONS']:
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
if 'preproc' in options['ACTIONS']:
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
else:
  yearMeanFiles = []
  for year in LOG['years']:
    yearMeanFiles.append(LOG['meansOfYears'][year])
dumpLog()
dbg(LOG)
# }}} ===================================================================================
# COMPUTE INITIAL VALUES FILES FOR LATER BIASES {{{ =====================================
LOG['init'] = cdo.seltimestep(1,input =  iFiles[0], output = '/'.join([options['ARCHDIR'],'%s_init.nc'%(options['EXP'])]))
# }}} ===================================================================================
# COMPUTE CELL MASK FOR LATER APPLICATION {{{ ===========================================
LOG['mask'] = cdo.selname('wet_c',input = '-seltimestep,1 %s'%(iFiles[0]), output = '/'.join([options['ARCHDIR'],'%s_mask.nc'%(options['EXP'])]))
# }}} ===================================================================================
# COMPUTE NUMBER OF VERTICAL LEVELS {{{ =================================================
LOG['depths'] = cdo.showlevel(input=LOG['mask'])[0].split()
# }}} ===================================================================================
# COMPUTE SINGLE YEARMEAN FILES {{{ =====================================================
ymFile = yearMeanFileName(options['ARCHDIR'],options['EXP'],LOG['years'][0],LOG['years'][-1])
if ( not os.path.exists(ymFile) or options['FORCE'] or hasNewFiles):
# if (os.path.exists(ymFile)):
#   os.remove(ymFile)
  cdo.cat(input=" ".join(yearMeanFiles),output=ymFile)
  # rm ymFiles
  #map(lambda x: os.remove(x),ymFiles)
#else:
# print("Use existing ymFile: "+ymFile)
# }}} ===================================================================================
# COMPUTE SINGLE MEAN FILE FROM THE LAST COMPLETE 30 YEARS {{{ =====================================================
lastYearsPeriod    = min(30,len(LOG['years'])/2)
lastYearsStartYear = -1*lastYearsPeriod - 2
lastYearsEndYear   = lastYearsStartYear + lastYearsPeriod

lastYearsFiles = []
for y in LOG['years'][lastYearsStartYear:lastYearsEndYear]:
    lastYearsFiles.append(LOG[str(y)])
yearMonMeanFile    = yearMonMeanFileName(options['ARCHDIR'],
                                         options['EXP'],
                                         LOG['years'][lastYearsStartYear],
                                         LOG['years'][lastYearsEndYear])
LOG['last30YearsMonMean']   = cdo.ymonmean(input = cdo.cat(input = ' '.join(lastYearsFiles),
                                                           output = yearMonMeanFile+'tmp'),
                                           output = yearMonMeanFile)
LOG['last30YearsMean']     = cdo.timmean(input = '-selyear,%s %s'%(','.join(LOG['years'][lastYearsStartYear:lastYearsEndYear]),ymFile),
                                         output = '%s/last30YearsMean_%s_%s-%s.nc'%(options['ARCHDIR'],
                                                                                    options['EXP'],
                                                                                    LOG['years'][lastYearsStartYear],
                                                                                    LOG['years'][lastYearsEndYear]))
LOG['last30YearsMeanBias'] = cdo.sub(input = ' %s %s'%(LOG['last30YearsMean'],LOG['init']),
                                     output = '%s/last30YearsMeanBias_%s_%s-%s.nc'%(options['ARCHDIR'],
                                                                                    options['EXP'],
                                                                                    LOG['years'][lastYearsStartYear],
                                                                                    LOG['years'][lastYearsEndYear]))
# }}} ===================================================================================
# PREPARE INPUT FOR PSI CALC {{{
# collect the last 20 years if there are more than 40 years, last 10 otherwise
if 'procPsi' in options['ACTIONS'] or 'plotPsi' in options['ACTIONS']:
  if len(LOG['years']) > 40:
    nyears4psi = 20
  else:
    nyears4psi = 10
  years4Psi = LOG['years'][-(nyears4psi+2):-1]
  dbg(LOG['years'])
  dbg(years4Psi)
  yearInfo  = '-'.join([years4Psi[0],years4Psi[-1]])
  psiModelVariableName = "u_vint_acc"
  psiGlobalFile = '/'.join([options['ARCHDIR'],'_'.join([psiModelVariableName,yearInfo])+'.nc'])
  cdo.timmean(input = "-selname,%s -selyear,%s/%s %s"%(psiModelVariableName,years4Psi[0],years4Psi[-1],ymFile),
              output = psiGlobalFile)
  psiSelectionConfig = {
          'indonesian_throughflow' : { 'lonlatbox' : '90,150,-20,40',},
          'gibraltar'              : { 'lonlatbox' : '-20,10,25,50',},
          'north_atlantic'         : { 'lonlatbox' : '-60,20,50,80',},
          'drake_passage'          : { 'lonlatbox' : '-90,-30,-80,-40',},
          'beringStrait'           : { 'lonlatbox' : '-180,-100,30,80',},
          'agulhas'                : { 'lonlatbox' : '10,50,-55,-15',},
          }

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
  regioCodes    = {'ArcticOcean' : 2, 'NorthAtlantic' : 4,'TropicalAtlantic' : 5, 'SouthernOcean' : 6}
  if 20 == len(LOG['depths']):
    regioDepths   = [110,215,895,2200]
  elif 40 == len(LOG['depths']):
    regioDepths   = [100,220,960,2290]

  else:
    regioDepths = [LOG['depths'][0],
                   LOG['depths'][len(LOG['depths'])/4],
                   LOG['depths'][len(LOG['depths'])/2],
                   LOG['depths'][len(LOG['depths'])-1]]
  regioVars     = ['t_acc','s_acc']
  regioMaskVar  = 'regio_c'
  dbg(regioDepths)
  dbg(LOG['depths'])
  dbg(len(LOG['depths']))
  regioMeanData = processRegionMean(options,LOG['mask'],LOG['init'],ymFile,regioCodes,regioDepths,regioVars,regioMaskVar)
# }}} ----------------------------------------------------------------------------------
# PREPARE INPUT FOR MOC PLOT {{{
if 'procMoc' in options['ACTIONS']:
  # collect all MOC files
  dbg(options['MOCPATTERN'])
  mocFiles        = sorted(glob.glob(options['MOCPATTERN']),key = mtime)
  mocFiles.pop(0)
  # skip the latest one if job is still running
  if options['JOBISRUNNING']:
    mocFiles.pop()
  # filter out emtpy files - this could happen when model finished before output timestep
  _mocFiles = []
  for mocfile in mocFiles:
    if (0 < os.stat(mocfile).st_size and 'Z' == mocfile[-1]):
      _mocFiles.append(mocfile)
  mocFiles = _mocFiles
  dbg(mocFiles)
  # default is to take the mean value ofthe at 10 years as input for the plotscript
  # this means 120 months, with monthly output, this is 120 timesteps
  mocNeededNSteps = 120
  mocLog          = {}
  scanFilesForTheirNTime(mocFiles,options['PROCS'],mocLog)
  dbg(mocLog)
  # check for the numbe rof timesteps in the last moc file
  mocLastNtime    = int(mocLog[mocFiles[-1]]) - 1 # avoid the last one, might be corrupted
  mocMeanFile     = '/'.join([options['ARCHDIR'],'mocMean'])
  if ( os.path.exists(mocMeanFile) ):
      os.remove(mocMeanFile)
  if mocNeededNSteps <= mocLastNtime:
    # take the last 120 values for timmeaninput
    mocMeanFile = cdo.timmean(input = "-seltimestep,%s/%s %s"%(mocLastNtime-mocNeededNSteps+1,mocLastNtime,mocFiles[-1]),
                              output = mocMeanFile)
  else:
    mocMeanFile = cdo.timmean(input = mocFiles[-1], output = mocMeanFile)
  dbg(mocMeanFile)
# }}} -----------------------------------------------------------------------------------
# PREPARE DATA FOR T,S,RHOPOT BIAS TO INIT {{{ ------------------------------------------
# target is a year mean file of fldmean data, but mean value computation
# should come at the very end of the processing chain
if 'procTSR' in options['ACTIONS']:
  t_s_rho_Input_1D = []
  t_s_rho_Input_2D = []
  for year in LOG['years']:
    t_s_rho_Input_1D.append(LOG[year])
  t_s_rho_Output_1D = computeMaskedTSRMeans1D(t_s_rho_Input_1D,['t_acc','s_acc','rhopot_acc'],
                                              LOG['init'],LOG['mask'],options['EXP'],options['ARCHDIR'],options['PROCS'])
  dbg(t_s_rho_Output_1D)
  t_s_rho_Output_2D = applyMask(' -selname,t_acc,s_acc,rhopot_acc %s'%(LOG['last30YearsMeanBias']),
                                LOG['mask'],
                                '/'.join([options['ARCHDIR'],'TSR_2D_%s_complete_timmean_masked.nc'%(options['EXP'])]))
  dbg(t_s_rho_Output_2D)
# }}} -----------------------------------------------------------------------------------
# DIAGNOSTICS ===========================================================================
# PSI {{{
if 'plotPsi' in options['ACTIONS']:
  plotFile = options['PLOTDIR']+'/'+"_".join(["psi",yearInfo,options['EXP'],options['TAG']+'.png'])
  title    = "Bar. Streamfunction for %s\n (file: %s)"%(options['EXP'],psiGlobalFile)
  cmd = '%s %s %s'%(options['CALCPSI'], psiGlobalFile, " DEBUG=1 WRITEPSI=true AREA=%s TITLE='%s' PLOT=%s"%(options['GRID'],title,plotFile))
  dbg(cmd)
  plotCommands.append(cmd)
  if subprocess.check_call(cmd,shell=True,env=os.environ):
    print("ERROR: CALCPSI failed")
  # plot special areas
  for area, selection in psiSelectionConfig.iteritems():
    title = "Selected Stream function for %s (%s)"%(area,options['EXP'])
    plotFile = options['PLOTDIR']+'/'+"_".join(["psi",area,options['EXP'],options['TAG']+'.png'])
    cmd = '%s %s %s'%(options['CALCPSI'], psiGlobalFile, " AREA=local TITLE='%s' PLOT=%s BOX=%s "%(title,plotFile,selection['lonlatbox']))
    dbg(cmd)
    plotCommands.append(cmd)
    if subprocess.check_call(cmd,shell=True,env=os.environ):
      print("ERROR: CALCPSI failed")
# }}} ----------------------------------------------------------------------------------
# HORIZONTAL PLOTS: t,s,u,v,abs(velocity) {{{
if 'plotHorz' in options['ACTIONS']:
  # A) last year mean
  horizontalConfig = {
    'varNames'      : ['t_acc','s_acc','h_acc','u_acc','v_acc'],
    'iFile'         : LOG['meansOfYears'][LOG['years'][-2]],
    'availableVars' : cdo.showname(input = LOG[LOG['years'][-2]])[0].split(' '),
    'sizeOpt'       : '-xsize=1200 -ysize=800',
    'title'         : '%s: last year mean '%(options['EXP']),
    'tag'           : 'lastYearMean',
    'limits'        : PlotConfig,
  }
  plotHorizontal(horizontalConfig,options,hasNewFiles)
  # B) last 30-year mean
  horizontalConfig = {
    'varNames'      : ['t_acc','s_acc','h_acc','u_acc','v_acc'],
    'iFile'         : LOG['last30YearsMean'],
    'availableVars' : cdo.showname(input = LOG['last30YearsMean'])[0].split(' '),
    'sizeOpt'       : '-xsize=1200 -ysize=800',
    'title'         : '%s: last 30-year-mean '%(options['EXP']),
    'tag'           : 'last30YearMean',
    'limits'        : PlotConfig,
  }
  plotHorizontal(horizontalConfig,options,hasNewFiles)
# }}} ----------------------------------------------------------------------------------
# THROUGH FLOWS / ONLINE DIAGNOSTICS {{{
# for global grid only
if ( 'global' == options['GRID'] ):
  diagnosticFiles = sorted(glob.glob(os.path.sep.join([LOG['dataDir'],"oce_diagnostics-*txt"])),key=mtime)
  if options['JOBISRUNNING']:
    diagnosticFiles.pop()
  diagnosticTable = plotOnlineDiagnostics(diagnosticFiles,options)
else:
  diagnosticTable = ''
# }}} ----------------------------------------------------------------------------------
# ATLANTIC X-Section: t,s,rhopot  {{{ ================================
# for global grid only
XSectionPlotConfig = PlotConfigBias
XSectionPlotConfig['s_acc']['minVar'] = '-1.0'
XSectionPlotConfig['s_acc']['maxVar'] = '1.0'
XSectionPlotConfig['t_acc']['minVar'] = '-5.0'
XSectionPlotConfig['t_acc']['maxVar'] = '5.0'
if 'plotX' in options['ACTIONS']:
  if ( 'global' == options['GRID'] ):
    for varname in ['t_acc','s_acc','rhopot_acc']:
      # plot last yearmean state
      oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_'+'_'.join([options['EXP'],options['TAG']])])
      iFile4XSection = LOG['last30YearsMean'] # take last 30 year mean file
      # mask by devision
      iFile4XSection = cdo.div(input  = '-selname,%s %s %s'%(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                               output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
      if ( 'rhopot_acc' == varname ):
        # substract 1000
        iFile4XSection  = cdo.subc(1000.0,
                                   input = '-selname,%s %s'%(varname,iFile4XSection),
                                   output = os.path.splitext(iFile4XSection)[0]+'_%s_subc1000.nc'%(varname))
      title = '%s: last 30 year mean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],
             '-resolution=r360x180',
             '-selPoints=150',
             '-rStrg="-"',
             '-makeYLinear',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      cmd = ' '.join(cmd)
      dbg(cmd)
      plotCommands.append(cmd)
      subprocess.check_call(cmd,shell=True,env=os.environ)
  
      # plot bias to initialization
      oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_BiasToInit'+'_'.join([options['EXP'],options['TAG']])])
      iFile4XSection = LOG['last30YearsMeanBias'] # take last 30 yearmean bias to init
      # mask by devision
      iFile4XSection = cdo.div(input  = '-selname,%s %s %s'%(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                               output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
      title = '%s: last 30 year mean bias to init '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],
             '-resolution=r360x180',
             '-selPoints=150',
             '-rStrg="-"',
             '-makeYLinear',
             '-minVar=%s'%(XSectionPlotConfig[varname]['minVar']),
             '-maxVar=%s'%(XSectionPlotConfig[varname]['maxVar']),
             '-numLevs=%s'%(XSectionPlotConfig[varname]['numLevs']),
             '-colormap=BlueWhiteOrangeRed',  # put white in the middle
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      cmd = ' '.join(cmd)
      dbg(cmd)
      plotCommands.append(cmd)
      subprocess.check_call(cmd,shell=True,env=os.environ)

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
#           '-oType='+options['OFORMAT'],
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
        ofile = '/'.join([options['PLOTDIR'],'_'.join(['.regioMean',location,varname,str(depth)+'m',options['EXP'],options['TAG']])+'.png'])

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
if 'plotMoc' in options['ACTIONS']:
  mocPlotCmd = '%s %s'%(options['MOCPLOTTER'],mocMeanFile)
  mocPlotSetup = {
    'TITLE' : 'ICON: %s : %s %s'%(options['EXP'],mocFiles[-1],'[last 10y mean]'),
    'IFILE' : mocMeanFile,
    'OFILE' : '/'.join([options['PLOTDIR'],'MOC_%s'%(options['EXP'])]),
    'OTYPE' : 'png',
  }
  for k,v in mocPlotSetup.iteritems():
    os.environ[k] = v
  plotCommands.append(mocPlotCmd)
  if subprocess.check_call(mocPlotCmd,shell=True,env=os.environ):
    print("ERROR: MOCPLOT failed")
  # environment cleanup
  for k in mocPlotSetup.keys():
    os.environ.pop(k)
# }}} ----------------------------------------------------------------------------------
# MIXED LAYER PLOT PLOT {{{
# take march and septempber from the last couple of years and just plot is
if ('procMld' in options['ACTIONS']):
  mldLastYear = LOG[LOG['years'][-2]]
  mldData     = cdo.timmax(input = ' -selname,mld,wet_c %s'%(mldLastYear), output = '/'.join([options['ARCHDIR'],'_'.join(['mld',os.path.basename(mldLastYear)])]))
if ('plotMld' in options['ACTIONS']):
  mldPlot = '/'.join([options['PLOTDIR'],'_'.join(['mld',os.path.basename(mldLastYear)])])
  mldCmd  = [options['ICONPLOT'],
      '-iFile=%s'%(mldLastYear),
      '-varName=mld',
      '-isIcon',
      '-colormap=WhiteBlueGreenYellowRed',
      '-maxVar=5000','-minVar=500',
      '-maskName=wet_c',
      '-oType='+options['OFORMAT'],'-tStrg="%s: ~C~ maximum MLD in the last year"'%(options['EXP']),
      '-oFile=%s'%(mldPlot)]
  print(mldCmd)
  mldCmd = ' '.join(mldCmd)
  if subprocess.check_call(mldCmd,shell=True,env=os.environ):
    print('CMD: %s has failed!'%(mldCmd))
  plotCommands.append(mldCmd)

  # nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$fslo -varName=mld -oType=png -colormap=WhiteBlueGreenYellowRed -timeStep=$t -oFile=mld_LIN_$t -maxVar=5000 -minVar=0 -maskName=wet_c
  # nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$fslo -varName=mld -oType=png -colormap=WhiteBlueGreenYellowRed -timeStep=$t -oFile=mld_LOG_$t -maxVar=5000 -minVar=0 -selMode=halflog -maskName=wet_c
# }}} ----------------------------------------------------------------------------------
# SEA ICE PLOT {{{
# take march and septempber from the last 30 year mean
if ('procIce' in options['ACTIONS']):
  iceData = cdo.selname('hi_acc,conc_acc',
                        input = LOG['last30YearsMonMean'],
                        output = '/'.join([options['ARCHDIR'],'sea_ice_%s'%(os.path.basename(LOG['last30YearsMean']))]))
  meanIceThicknessData = cdo.setunit('m',
                                     input = " -expr,'hiMean=hi_acc*conc_acc' %s"%(iceData),
                                     output = '/'.join([options['ARCHDIR'],'meanIceThickness_%s'%(os.path.basename(iceData))]))
if ('plotIce' in options['ACTIONS']):
  for mapTypeIndex,mapType in enumerate(['NHps','SHps']):
    month = [3,9][mapTypeIndex]
    iFile = cdo.selmon(month,
                       input = iceData,
                       output = '/'.join([options['ARCHDIR'],'month%s_%s'%(month,os.path.basename(iceData))]))
    icePlotConfig = {
      'varNames'      : ['hi_acc','conc_acc'],
      'iFile'         : iFile,
      'availableVars' : cdo.showname(input = iceData)[0].split(' '),
      'sizeOpt'       : '-xsize=1200 -ysize=800',
      'title'         : '%s:~C~ SEA ICE month %s, last-30-yearmonmean'%(options['EXP'],month),
      'tag'           : 'last30YearMeanICE_%s'%(mapType),
      'mapType'       : mapType,
      'colormap'      : 'WhBlGrYeRe',
      'maskFile'      : LOG['mask'],
    }
    plotHorizontal(icePlotConfig,options,hasNewFiles)
    iFile = cdo.selmon(month,
                       input = meanIceThicknessData,
                       output = '/'.join([options['ARCHDIR'],'month%s_%s'%(month,os.path.basename(meanIceThicknessData))]))
    icePlotConfig = {
      'varNames'      : ['hiMean'],
      'iFile'         : iFile,
      'availableVars' : cdo.showname(input = meanIceThicknessData)[0].split(' '),
      'sizeOpt'       : '-xsize=1200 -ysize=800',
      'title'         : '%s:~C~ SEA ICE, month %s, last-30-yearmonmean ~C~ mean sea ice thickness'%(options['EXP'],month),
      'tag'           : 'last30YearMeanICE_%s'%(mapType),
      'mapType'       : mapType,
      'colormap'      : 'WhBlGrYeRe',
      'maskFile'      : LOG['mask'],
    }
    plotHorizontal(icePlotConfig,options,hasNewFiles)

# }}} ----------------------------------------------------------------------------------
# take march and septempber from the last couple of years and just plot is
# T S RHOPOT BIAS PLOT {{{
if 'plotTSR' in options['ACTIONS']:
  # global mean bias over depth and time
  t_s_rho_PlotSetup = PlotConfigBias
  for varname in t_s_rho_PlotSetup.keys():
    oFile = '/'.join([options['PLOTDIR'],varname+'_biasToInit_inDepth_overTime'+'_'.join([options['EXP'],options['TAG']])])
    title = '%s: %s bias to init '%(options['EXP'],varname)
    cmd = [options['ICONPLOT'],
           '-iFile=%s'%(t_s_rho_Output_1D),
           '-varName=%s'%(varname),
           '-oType='+options['OFORMAT'],
           '-hov',
           '-rStrg="-"',
           '-bStrg="-"',
           '-lStrg=%s'%(varname),
  #        '-withLineLabels',
           '+withLines',
           '-makeYLinear',
           '-minVar=%s'%(t_s_rho_PlotSetup[varname]['minVar']),
           '-maxVar=%s'%(t_s_rho_PlotSetup[varname]['maxVar']),
           '-numLevs=%s'%(t_s_rho_PlotSetup[varname]['numLevs']),
           '-colormap=BlWhRe',
    #      '-tStrg="%s"'%(title),
           '-oFile=%s'%(oFile)]
    cmd = ' '.join(cmd)
    dbg(cmd)
    plotCommands.append(cmd)
    if subprocess.check_call(cmd,shell=True,env=os.environ):
      print('CMD: %s has failed!')

# }}} ----------------------------------------------------------------------------------
# FINAL DOCUMENT CREATION {{{ ===========================================================
if 'finalDoc' in options['ACTIONS']:
  createOutputDocument(options['PLOTDIR'],diagnosticTable,'_'.join(['ALL',options['EXP'],options['TAG']]),options['DOCTYPE'],options['DEBUG'])
# }}} ----------------------------------------------------------------------------------
for plotcmd in plotCommands:
  print(plotcmd)

doExit
#
#
# vim:fdm=marker
# vim:sw=2
