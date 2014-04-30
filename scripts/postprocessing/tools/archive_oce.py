#!/usr/bin/env python
import os,sys,glob,shutil,json,subprocess,multiprocessing
from cdo import *
import matplotlib
matplotlib.use('Agg')
import pylab,numpy,dateutil.parser


# ==============================================================================
# OPTION HANDLING {{{ =========================================================
options = {
            'ARCHDIR'     : './archive',
            'EXP'         : 'oce_mpiom',
            'FILEPATTERN' : 'oce_mpiom/oce_mpiom_*.nc*',
            'DEBUG'       : False,
            'FORCE'       : False,
            'CALCPSI'     : '/scratch/mpi/CC/mh0287/users/m300064/src/icon-HEAD/scripts/postprocessing/tools/calc_psi.py',
            'TAG'         : 'r1xxxx',
            'ICONPLOT'    : 'nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools',
            'PROCS'       : 8,
           }

optsGiven = sys.argv[1:]
for optVal in optsGiven:
    key,value    = optVal.split('=')
    if key in ['FORCE','DEBUG']:
        value = value.lower() in ['true','1']
    if 'PROCS' == key:
      value = int(value)

    options[key] = value
# }}}
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

"" " load the last archiving status ir available """
def loadLog():
  if os.path.exists(LOGFILE):
    dbg("Load last status from LOGFILE")
    LOG = json.load(open(LOGFILE))
  else:
    dbg("Could not find file logile: '#{@logile}'! Create a new one ....")
    options['FORCE'] = True
    LOG = {}

  return LOG

""" get year from yearly output file """
def yearFrom(yfilename):
  return os.path.basename(yfilename).split('.')[0].split('_')[-1]

""" append list of image together """
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
  map(lambda x: os.remove(x),imageMap.values())

""" warpper around cdo showyear for multiprocessing """
def showyear(file):
  return cdo.showyear(input = file)[0].split(' ')

""" collect the years, which aree stored in to files """
def scanFilesForTheirYears(fileList,procs,log):

  log['yearsOfFiles'] = {}
  pool = multiprocessing.Pool(14)
  lock = multiprocessing.Lock()
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
# }}}
# INTERNALS {{{ ================================================================
LOGFILE         = 'archive.log'
LOG             = loadLog();dbg(LOG)
cdo             = Cdo()
cdo.debug       = options['DEBUG']
cdo.forceOutput = options['FORCE']
# }}}
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
def usage():
  return """
  # USAGE =====================================================================
  #
  #   ./archive.py EXP=<expname> FILEPATTERN=<pattern> ARCHDIR=<dir> TAG=<rxxxxx>
  #
  # ===========================================================================
  #
  # expname     : string for naming things
  # tag         : additional placeholder, e.g. for revisions, if not used, set it to ''
  # filepattern : quoted wildcard for model result relative to the 'experiments' dir
  # archdir     : directory, where the archived data is placed (default: ./archive)
  #
  # more options are (keys area llways uppercase)
  # DEBUG       : False/True
  # CALCPSI     : path, where to find the psi calculation and plotting scipt
  # ICONPLOT    : call for icon_plot.ncl, default is 'nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools'
  #
  # #
  # !!!! THIS SCRIPT IS MEANT TO BE STARTED IN THE MODEL RESULTS DIRECTORY !!!!
  """
#=============================================================================== 
# Check input:
#   options
dbg(options)
#   are the files accessable?
ifiles = glob.glob(options['FILEPATTERN'])
dbg(ifiles); LOG['ifiles'] = ifiles
if 0 == len(ifiles):
  print usage()
  print("Could not find any result files!")
  exit(1)
else:
  print("will use \n" + "\n".join(ifiles))
# is the archive dir accessable
if not os.path.isdir(options['ARCHDIR']):
  os.makedirs(options['ARCHDIR'])
# =======================================================================================
# DATA SPLITTING {{{ ====================================================================
#TODO:
# scanFilesForTheirYears(ifiles,options['PROCS'],LOG)
# dbg(LOG)
# doExit()

if (options['FORCE'] or ('splityear?' in LOG.keys() and (False == LOG['splityear?']))):
  try:
       
    if os.path.isdir(options['ARCHDIR']):
      dbg("remove ARCHDIR: "+options['ARCHDIR'])
      shutil.rmtree(options['ARCHDIR'])
      os.makedirs(options['ARCHDIR'])

    cdo.splityear(input  = "-cat '%s'"%(options['FILEPATTERN']),
                  output = '/'.join([options['ARCHDIR'],options['EXP']+'_']),
                  options = '-v')
    LOG['splityear?'] = True
    LOG['splitfiles'] = glob.glob('/'.join([options['ARCHDIR'],options['EXP']+'_'])+'*')

    # some error handling: cdo puts a '*' at the end, if its at the end of the pattern'
    for i,f in enumerate(LOG['splitfiles']):
      if '*' == f[-1]:
        shutil.move(f,f[0:-1])
        LOG['splitfiles'][i] = f[0:-1]
    LOG['years'] = map(lambda x: yearFrom(x), LOG['splitfiles'])
    for year in LOG['years']:
      yearlyFile  = '/'.join([options['ARCHDIR'],options['EXP']+'_'])+year+'.nc'
      if os.path.exists(yearlyFile):
        LOG[year] = yearlyFile

    dumpLog()

  except:
    print("splityear failed somehow")
    LOG['splityear?'] = False
    doExit(1)

dbg(LOG)
# }}} ===================================================================================
# COMPUTE YEARMEAN FILES {{{ ============================================================
ymFile = '/'.join([options['ARCHDIR'],'_'.join([options['EXP'],'yearmean.nc'])])
if not os.path.exists(ymFile):
  def _createYearMeanOf(ifile,ofile):
    cdo.yearmean(input=ifile,output=ofile)

  ymFiles = []
  pool  = multiprocessing.Pool(options['PROCS'])
  for f in LOG['splitfiles']:
    ofile = '/'.join([options['ARCHDIR'],"ym_"+os.path.basename(f)])
    pool.apply_async(_createYearMeanOf,[f,ofile])
    ymFiles.append(ofile)

  pool.close()
  pool.join()

  cdo.cat(input=" ".join(ymFiles),output=ymFile)
  # rm ymFiles
  map(lambda x: os.remove(x),ymFiles)
else:
  print("Use existing ymFile: "+ymFile)
# }}} ===================================================================================
# PREPARE INPUT FOR PSI CALC {{{
# collect the last 20 years if there are more than 40 years, last 10 otherwise
if len(LOG['years']) > 40:
  nyears4psi = 20
else:
  nyears4psi = 10
years4Psi = LOG['years'][-(nyears4psi+2):-1]
dbg(years4Psi)
uvintName = 'u_vint_acc'
yearInfo = '-'.join([years4Psi[0],years4Psi[-1]])
uvintFile = '/'.join([options['ARCHDIR'],'_'.join([uvintName,yearInfo])+'.nc'])
cdo.timmean(input = "-selname,%s -selyear,%s/%s %s"%(uvintName,years4Psi[0],years4Psi[-1],ymFile),
            output = uvintFile)
# }}}
# PREPARE INPUT FOR PROFILES {{{
varNames = ['t_acc','s_acc','u_acc','v_acc']
#  collect the data first (first 10 years)
varFiles = []
for year in LOG['years'][0:10]:
# varFiles.append(cdo.selname(' '.join(varNames),
#   input=" -remapnn,lon=-30.0_lat=-65.0 %s"%(LOG[year]),
#   output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))
  varFiles.append(cdo.remapnn('lon=-30.0_lat=-65.0' , 
                              input = ' -selname,%s %s'%(','.join(varNames),LOG[year]),
                              output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))

varData = cdo.cat(input=' '.join(varFiles),
                    output="%s/_data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP']))
# add absolute velocity
velocityData = cdo.expr("'velocity=sqrt(u_acc*u_acc+v_acc*v_acc);'",input=varData,output='/'.join([options['ARCHDIR'],
                                                                                                   'vel4profile.nc']))
profileData = cdo.merge(input='%s %s'%(varData,velocityData),output="%s/data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP']))
# }}}
# PREPARE INPUT FOR REGIONAL MEANS from yearMean output {{{
# helper for parallelization
def _computeRegioMean(depth,varname,ifile,vertMask,regioMask,ofile):
  # mask out region
  # vertical interpolation to target depth
  # mean value computaion
  cdo.fldmean(input = '-mul -div -sellevel,%s -selname,%s %s %s %s'%(depth,varname,ifile,vertMask,regioMask), output = ofile)

# setup
regioCodes    = {'NorthAtlantic' : 4,'TropicalAtlantic' : 5, 'SouthernOcean' : 6}
regioDepths   = [110,215,895,2200]
regioVars     = ['t_acc','s_acc']
regioMaskVar  = 'regio_c'
regioMeanData = {}
regioPool     = multiprocessing.Pool(options['PROCS'])
regioLock     = multiprocessing.Lock()

# create the regio mask first
regioMasks    = {}
for location, regioCode in regioCodes.iteritems():
  ofile     = '/'.join([options['ARCHDIR'],'_'.join(['regioMask',location])+'.nc'])
  ofileTemp = '/'.join([options['ARCHDIR'],'_'.join(['_regioMask',location])+'.nc'])
  cdo.eqc(regioCode,input = '-selname,%s -seltimestep,1 %s'%(regioMaskVar,LOG['ifiles'][0]),output = ofileTemp)
  regioMasks[location] = cdo.div(input = '%s %s'%(ofileTemp,ofileTemp),output = ofile)
# create the mask from 3d mask wet_c
regioVertMasks    = {}
for depth in regioDepths:
  depth = str(depth)
  ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioVertMask',depth+'m'])+'.nc'])
  regioVertMasks[depth] = cdo.sellevel(depth,input = '-div -selname,wet_c -seltimestep,1 %s -selname,wet_c -seltimestep,1 %s'%(LOG['ifiles'][0],LOG['ifiles'][0]),output = ofile)
# compute the regional mean values
for location, regioCode in regioCodes.iteritems():
  regioMeanData[location] = {}
  for depth in regioDepths:
    regioMeanData[location][str(depth)] = {}
    for varname in regioVars:
      regioMeanData[location][str(depth)][varname] = {}
      ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,str(depth)+'m'])+'.nc'])
      regioPool.apply_async(_computeRegioMean,[depth,varname,ymFile,regioVertMasks[str(depth)],regioMasks[location],ofile])
      regioLock.acquire()
      regioMeanData[location][str(depth)][varname] = ofile 
      regioLock.release()
regioPool.close()
regioPool.join()
# }}}
# DIAGNOSTICS ===========================================================================
# PSI {{{
plotFile = options['ARCHDIR']+'/'+"_".join(["psi",yearInfo,options['EXP'],options['TAG']+'.png'])
if not os.path.exists(plotFile):
  cmd = '%s %s %s'%(options['CALCPSI'], uvintFile, "LEVELS=20 PLOT="+plotFile)
  dbg(cmd)
  if subprocess.check_call(cmd,shell=True,env=os.environ):
    print("ERROR: CALCPSI failed")
# }}}
# DRAKE FLOW {{{
# cat all diagnostics together
diagnosticFiles  = glob.glob("oce_diagnostics-*txt")
diagnosticJoined = options['ARCHDIR']+'/'+'_'.join([options['EXP'],"diagnostics.txt"])

header     = open(diagnosticFiles[0],'r').readline().split(' ')
dateIndex, drakeIndex  = header.index('date'), header.index('drake_passage')

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
dates = all_diag[dateIndex]
drake = all_diag[drakeIndex]

values = map(float,numpy.array(drake)[1:-1])
dates  = map(dateutil.parser.parse,dates[1:-1])
# complete timeseries
ofile = "%s/drake_complete_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG'])
if ( not os.path.exists(ofile) or options['FORCE']):
  pylab.title("%s , %s :Drake TF - all %s years"%(options['EXP'],options['TAG'],len(LOG['years'])),fontsize=8)
  pylab.grid()
  pylab.plot_date(dates, values, linestyle='-',marker='.')  
  pylab.savefig(ofile)
  pylab.clf() 
# first 10 years
ofile = "%s/drake_first10Years_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG'])
if ( not os.path.exists(ofile) or options['FORCE']):
  pylab.title("%s , %s :Drake TF - first 10y"%(options['EXP'],options['TAG']),fontsize=8)
  pylab.grid()
  pylab.plot_date(dates[0:120], values[0:120], linestyle='-',marker='.')  
  pylab.savefig(ofile)
  pylab.clf() 
# first 20 years
ofile = "%s/drake_first20Years_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG'])
if ( not os.path.exists(ofile) or options['FORCE']):
  pylab.title("%s , %s :Drake TF - first 20y"%(options['EXP'],options['TAG']),fontsize=8)
  pylab.grid()
  pylab.plot_date(dates[0:240], values[0:240], linestyle='-',marker='.')  
  pylab.savefig(ofile)
# }}}
# SOUTH OCEAN t,s,y,v profile at 30w, 65s  {{{ ================================
#  create hovmoeller-like plots
for varname in ['t_acc','s_acc','u_acc','v_acc','velocity']:
  # run icon_plot.ncl
  oFile = '/'.join([options['ARCHDIR'],varname+'_profile_30w-65s_'+'_'.join([options['EXP'],options['TAG']])])
  if ( not os.path.exists(oFile+'.png') or options['FORCE']):
    title = '%s: 10y profile at 30W,65S'%(options['EXP'])
    cmd = [options['ICONPLOT'],
           '-iFile=%s'%(profileData),
           '-hov',
           '-varName=%s'%(varname),
           '-oType=png',
           '-rStrg="-"',
           '-tStrg="%s"'%(title),
           '-oFile=%s'%(oFile)]
    dbg(' '.join(cmd))
    subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)
# }}}
# REGIO MEAN PROFILES {{{ ================================
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
      ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,str(depth)+'m',options['EXP'],options['TAG']])+'.png'])

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
                            '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,options['EXP'],options['TAG']+'.png'])])) 
# }}}
# vim:fdm=marker
