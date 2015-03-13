# ================================================================================
function diffWithExitCode {
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
  refDir=$1
  expDir=$2

  refList=`ls ${refDir}/*`
  for refFile in ${refList}; do 
    refFileBasename=$(basename ${refFile})
    case "${refFileBasename}" in
      *.nc*)
        DIFF='diffWithExitCode'
        ;;
      *)
        DIFF='diff'
        ;;
    esac
    ${DIFF} ${refFile} ${expDir}/${refFileBasename}
  done
}
