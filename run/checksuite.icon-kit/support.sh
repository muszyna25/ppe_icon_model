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
	 icon_data_poolFolder=
     ;;
   xfh2*) :
     echo "...FH2  at KIT"; CENTER="IMK"
     input_folder="/pfs/imk/ICON/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
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
complete_output_folder=${output_folder}/${d}/${exp_name}
OUTDIR=$complete_output_folder
OUTDIR_PREFIX=`dirname ${OUTDIR}`
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run

cat > $output_script << EOF
#!/bin/bash
CENTER=$CENTER
basedir=$ICON_FOLDER
icon_data_poolFolder=$icon_data_poolFolder
EXPNAME=atm_amip_test_kit
OUTDIR=$complete_output_folder
ICONFOLDER=$ICON_FOLDER
ARTFOLDER=$ART_FOLDER
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

ln -sf $ICON_FOLDER/data/rrtmg_lw.nc rrtmg_lw.nc
ln -sf $ICON_FOLDER/data/ECHAM6_CldOptProps.nc ECHAM6_CldOptProps.nc

ln -sf ${ARTFOLDER}/runctrl_examples/init_ctrl/mozart_coord.nc ${OUTDIR}/mozart_coord.nc
ln -sf ${ARTFOLDER}/runctrl_examples/init_ctrl/Linoz2004Br.dat ${OUTDIR}/Linoz2004Br.dat
ln -sf ${ARTFOLDER}/runctrl_examples/init_ctrl/Simnoy2002.dat ${OUTDIR}/Simnoy2002.dat

EOF

}
function create_footer
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/$1.run

cat >> $output_script << EOF
	
cp -p $ICON_FOLDER/build/x86_64-unknown-linux-gnu/bin/icon ./icon.exe
	
echo $CENTER
	case x"$CENTER" in
      xFZJ*)
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
;;
	esac

EOF
}
function create_body

{
 cat ${ART_FOLDER}/runctrl_examples/run_scripts//exp.testsuite.amip_base >> $output_script
}

