#!/usr/bin/env python
import os,sys,glob,shutil,json,subprocess,multiprocessing
from cdo import *

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
# DIAGNOSTICS ===========================================================================
# PSI {{{
plotFile = options['ARCHDIR']+'/'+"_".join(["psi",yearInfo,options['EXP'],options['TAG']+'.png'])
dbg(plotFile)
if not os.path.exists(plotFile):
  if os.system('%s %s %s'%(options['CALCPSI'], uvintFile, "LEVELS=20 PLOT="+plotFile)):
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

import pylab,numpy,dateutil.parser

values = map(float,numpy.array(drake)[1:-1])
dates  = map(dateutil.parser.parse,dates[1:-1])

pylab.title("%s , %s :Drake TF - all %s years"%(options['EXP'],options['TAG'],len(LOG['years'])))
pylab.grid()
pylab.plot_date(dates, values, linestyle='-',marker='.')  
pylab.savefig("%s/drake_complete_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG']))


pylab.clf() 
pylab.title("%s , %s :Drake TF - first 10y"%(options['EXP'],options['TAG']))
pylab.grid()
pylab.plot_date(dates[0:120], values[0:120], linestyle='-',marker='.')  
pylab.savefig("%s/drake_first10Years_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG']))

pylab.clf() 
pylab.title("%s , %s :Drake TF - first 20y"%(options['EXP'],options['TAG']))
pylab.grid()
pylab.plot_date(dates[0:240], values[0:240], linestyle='-',marker='.')  
pylab.savefig("%s/drake_first20Years_%s_%s.png"%(options['ARCHDIR'],options['EXP'],options['TAG']))

# }}}
# SOUTH OCEAN t,s,y,v profile at 30w, 65s  {{{ ================================
#  create hovmoeller-like plots
for varname in ['t_acc','s_acc','u_acc','v_acc','velocity']:
  # run icon_plot.ncl
  oFile = '/'.join([options['ARCHDIR'],varname+'_profile_30w-65s_'+'_'.join([options['EXP'],options['TAG']])])
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
  proc = subprocess.Popen(' '.join(cmd),
                              shell  = True,
                          stderr = subprocess.PIPE,
                          stdout = subprocess.PIPE,
                          env    = os.environ)
  retvals = proc.communicate()
  ret =  {"stdout"     : retvals[0]
          ,"stderr"     : retvals[1]
          ,"returncode" : proc.returncode}
  print(ret['stdout'])
  print(ret['stderr'])
# }}}
# vim:fdm=marker
