# ================================================================================
function accumulationDiff {
  ifile=$1
  varName=$2
  maskName=$3
  ofile=`mktemp`

  cdo -s diffv -div -selname,$varName $ifile -selname,$maskName $ifile -div -selname,${varName}_acc $ifile -selname,$maskName $ifile > $ofile

  nDiff=$(wc -l < $ofile)

  if (( $nDiff > 0 )); then
    cat $ofile;
  fi

  rm -f $ofile
  return $nDiff
}

function directoryDiff {
  refDir=$1
  expDir=$2

  for refFile in ${refDir}/*; do 
    refFileBasename=$(basename ${refFile})
    case "${refFileBasename}" in
      *.nc*)
        DIFF='cdo diffv'
        ;;
      *)
        DIFF='diff'
        ;;
    esac
    ${DIFF} ${refFile} ${expDir}/${refFileBasename}
  done
}
