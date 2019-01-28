#!/bin/ksh 

function set_cluster {
 ##
 echo " Running on " >&2
 HPC=`(hostname || uname -n) 2>/dev/null | sed 1q`
 case x"$HPC" in #(
   xuc1*) :
     echo "...UC1 at KIT"; CENTER="IMK" 
     input_folder="/pfs/imk/ICON/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
	 icon_data_poolFolder=/pfs/work6/workspace/scratch/ln1297-AMIP_release2.5-0/MISTRAL
     ;;
   xfh2*) :
     echo "...FH2  at KIT"; CENTER="IMK"
     input_folder="/pfs/imk/ICON/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
	 icon_data_poolFolder=/pfs/work6/workspace/scratch/ln1297-AMIP_release2.5-0/MISTRAL

     ;;
   xxce*) :
     echo "...XCE at DWD"; CENTER="DWD"
	 FILETYPE="4" 
	 if [[ ":${PE_ENV}" != ':CRAY' ]] ; then
          module unload libdwd grib_api
          module unload cray-netcdf netcdf perftools stat
          module unload craype-hugepages2M
          case "{PE_ENV}" in
             GNU)   prgenv_="PrgEnv-gnu"   ;;
             INTEL) prgenv_="PrgEnv-intel" ;;
          esac
          module swap ${prgenv_} PrgEnv-cray
       fi
       module swap cce cce/${compiler_version}
       module load cray-netcdf perftools stat
       module load libdwd
       module load grib_api
       module list
       ##       
       ;;
   xlce*) :
     echo "...LCE at DWD"; echo "  ERROR: Testsuite has to be run and compiled on XCE !!!", exit 201 ;; #(
   xjuwels*) :
     echo "...juwels at FZJ"; CENTER="FZJ" 
     icon_data_poolFolder=/gpfs/homea/hka21/hka211/INPUT/AMIP/atm_amip_base/
	 FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
	 ;; #(
   xmlogin*) :
     echo "...MISTRAL at DKRZ"; CENTER="DKRZ"
     icon_data_poolFolder=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/ 
     input_folder="~/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="/scratch/b/${USER}/TESTSUITE_OUTPUT/"
	 ;;
   *) :
     echo "...unknown HPC" ; exit 202 ;; #(
 esac
 ##
 ICON_FOLDER=$(pwd)/../..
 ART_FOLDER=$ICON_FOLDER/src/art


 }

function read_setup {
#declare -a experiment
#declare -a script
#awk -va="experiment=$experiment; script=$script" 'BEGIN{n=1;}
#     { FS=" "; $0=$0; script[n]=$1; experiment[n]=$2;
#		 print n $1 script[n];
#
#     }' 'setup_file.txt' 
#awk -F '{ print $1, $2}' 'setup_file.txt' | read var1 var2
#echo $var1 ' ' $var2
i=0

for line in `awk '{print $2}' $1`; do
  experiment[$i]=$line
  i=`expr $i + 1`
done
i=0

for line in `awk '{print $1}' $1`; do
  script[$i]=$line
  i=`expr $i + 1`
done

i=0
for line in `awk '{print $3}' $1`; do
  config[$i]=$line
  i=`expr $i + 1`
done

i=0
for line in `awk '{print $4}' $1`; do
  walltime[$i]=$line
  i=`expr $i + 1`
done


i=0
for line in `awk '{print $5}' $1`; do
  queue[$i]=$line
  i=`expr $i + 1`
done


i=0
for line in `awk '{print $6}' $1`; do
  nodes[$i]=$line
  i=`expr $i + 1`
done



number_of_experiments=$i

 }

 function usage
{
    echo "This is a script for creating buildbot / testsuite / experiment scripts"
    echo "for different clusters"
    echo "Usage:"
    echo "./create_testsuite_runs --input_file=YOUR_FILE"
    echo ""
}

function create_header
{
d=`date -d today +%Y%m%d`
complete_output_folder=${output_folder}/${d}/$1
OUTDIR=$complete_output_folder
OUTDIR_PREFIX=`dirname ${OUTDIR}`
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run
OUTDIR=$complete_output_folder
icon_data_poolFolder=$icon_data_poolFolder
EXPERIMENT=$1
lart=$lart

cat > $output_script << EOF
#!/bin/bash
CENTER=$CENTER
basedir=$ICON_FOLDER
icon_data_poolFolder=$icon_data_poolFolder
EXPNAME=atm_amip_test_kit
OUTDIR=$complete_output_folder
ICONFOLDER=$ICON_FOLDER
ARTFOLDER=$ART_FOLDER
INDIR=$input_folder
EXP=$EXPERIMENT
lart=$lart

FILETYPE=4
COMPILER=intel
restart=.False.
read_restart_namelists=.False.



# Remove folder ${EXP} from OUTDIR for postprocessing output
OUTDIR_PREFIX=`dirname ${OUTDIR}`

# Create output directory and go to this directory

if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

cd $OUTDIR


EOF

}
function create_footer
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run

cat >> $output_script << EOF
	
cp -p $ICON_FOLDER/build/x86_64-unknown-linux-gnu/bin/icon ./icon.exe
EOF
 case x"$HPC" in #(
      xjuwels*)
cat >> $output_script << EOF
	
cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --nodes=4
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=00:15:00
#SBATCH --partition=batch
module load Intel
module load ParaStationMPI
module load Python/2.7.14
module load netCDF-Fortran/4.4.4
module load netCDF/4.6.1
srun ./icon.exe

ENDFILE


chmod +x job_ICON
sbatch job_ICON

EOF
;;

   xfh2*)
cat >> $output_script << EOF
	   
cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --$3
#SBATCH --time=$2
#SBATCH --ntasks-per-node=20
#SBATCH --partition=$4

export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/pfs/imk/ICON/LIBRARIES_IFORT16/szip/lib:/pfs/imk/ICON/LIBRARIES_IFORT16/grib-api/lib

$5 

mpirun --bind-to core --map-by core --report-bindings ./icon.exe

ENDFILE

chmod +x job_ICON
sbatch job_ICON

EOF
;;
esac


chmod +x $output_script
}
function create_body

{
 cat ${ART_FOLDER}/$1 >> $output_script
}


function create_lart_loop_start
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run
cat >> $output_script << EOF
	
for iter in 1 2 ; do


if [ \$iter -eq 1 ] ; then

lart=.False.
if [ ! -d $OUTDIR/wo_ART ]; then
    mkdir -p $OUTDIR/wo_ART
fi

cd $OUTDIR/wo_ART
OUTDIR=$OUTDIR/wo_ART
else

lart=.True.
if [ ! -d $OUTDIR/w_ART ]; then
    mkdir -p $OUTDIR/w_ART
fi

cd $OUTDIR/w_ART

OUTDIR=$OUTDIR/w_ART

fi

EOF
}

function create_lart_loop_end
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run
cat >> $output_script << EOF

done

EOF
}

function check_action
{
output_script=$2.run
	
	case x"$1" in 
		xrun*)
          ./$output_script
		  ;;
		*none*)
		  ;;
  esac
}

function read_configure
{
 . ./${pwd}/../../config/set-up.info

output="module load ${use_load_modules}"
#i=0
#for element in $(echo ${use_load_modules}|sed -e " ");do
#   module_str[$i]="module load ${element}"
#  i=`expr $i + 1`
#done
#
#for item in ${!module_str[*]};do
#	 printf " %s \n" "${module_str[${item}]}"
#done

 }

