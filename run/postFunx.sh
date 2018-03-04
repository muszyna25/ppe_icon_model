# ================================================================================
function diffWithExitCode {
#   set -x 
  fileA=$1
  fileB=$2
  ofile=`mktemp`

  cdo -s diffv ${fileA} ${fileB} > $ofile

  nDiff=$(wc -l < $ofile)

  if (( $nDiff > 0 )); then
    cat $ofile;
  fi

  rm -f $ofile
  return $nDiff
}
function accumulationDiff {
  ifile=$1
  varName=$2
  maskName=$3

  diffWithExitCode "-div -selname,$varName $ifile -selname,$maskName $ifile " "-div -selname,${varName}_acc $ifile -selname,$maskName $ifile"
}

function directoryDiff {
#   set -x 
  refDir=$1
  expDir=$2
  nDiff=0

  if [[ ! -d $refDir ]]; then
    return 99
  fi


  refList=`ls ${refDir}/*`
  for refFile in ${refList}; do 
    if [[ -f $refFile ]]; then
      refFileBasename=$(basename ${refFile})
      case "${refFileBasename}" in
      *.nc*)
        DIFF='diffWithExitCode'
      ;;
      *)
        DIFF='true' # disabled for not nc-files
      ;;
      esac
      ${DIFF} ${refFile} ${expDir}/${refFileBasename}
      if (( $nDiff > 0 )); then
        return $nDiff
      fi
    fi
  done
}
