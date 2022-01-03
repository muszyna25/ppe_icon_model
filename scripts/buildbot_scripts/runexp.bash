#! /bin/bash
#------------------------------------------------------------------------------
# This script is called on all buildbot slaves and should be used to run the runscripts
#------------------------------------------------------------------------------

# Define of the text which is written in case of a script fail or OK run

TEXT_RUN_FAILED="FAILED"
TEXT_RUN_OK="OK"

#==============================================================================
stop_on_error()
{
    echo STATUS_IN_FILE ${STATUS_IN_FILE}

    # The first argument has priority if it is not 0:
    if test "x$1" != x0; then _status=$1; else _status=${STATUS_IN_FILE}; fi

    if test "x${_status}" = x0; then
      # WARNING: cover the case when ${_status} is an empty string or a
      # space-only string if you decide to refactor the condition above.
      :
    elif test 1 -le "$_status" 2>/dev/null && test 255 -ge "$_status" 2>/dev/null; then
      # ${_status} is a valid exit code in range [1;255]:
      printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
      exit ${_status}
    else
      # We avoid exiting with an invalid exit code (e.g. 256 or '')
      printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
      exit 255
    fi

    printf '%-50s : %s\n' $2 ${TEXT_RUN_OK} >> ../${LOOP_STATUS_FILE}
}

#------------------------------------------------------------------------------------

warning_on_error()
{
    # The first argument has priority if it is not 0:
    if test "x$1" != x0; then _status=$1; else _status=${STATUS_IN_FILE}; fi

    if test "x${_status}" = x0; then
      # WARNING: cover the case when ${_status} is an empty string or a
      # space-only string if you decide to refactor the condition above.
      printf '%-50s : %s\n' $2 ${TEXT_RUN_OK} >> ../${LOOP_STATUS_FILE}
      if [ "${build_post_file}" == "true" ]
      then
        name="`echo $2 | cut -d '.' -f2`"
        RUN_POST="${RUN_POST} post.${name}.run"
        echo "${RUN_POST}" > ./run_post_list

        compare_run="post.${name}_compare_restarts.run"
        if [  -a ./${compare_run} ] ; then
          RUN_POST_COMP="${RUN_POST_COMP} $compare_run"
          echo "${RUN_POST_COMP}" > ./run_post_comp_list
        fi
      fi
    else
      printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
      echo "*********** WARNING: script failed  ${_status} ****************"
    fi
}

#------------------------------------------------------------------------------------

return_ok()
{
    # Exit with status = 0 = OK
    # Arguments:
    #   $1 =  message
    echo "return_ok()"
    echo "$1"

    exit 0
}

#------------------------------------------------------------------------------------



#==============================================================================
# runs the scripts 
run_scripts_submit()
{

  echo "|======================================================|"
  echo "|                                                      |"
  echo "|         Running exp-scripts                          |"
  echo "|                                                      |"
  echo "|======================================================|"

  LOOP_STATUS_FILE="LOOP_STATUS_EXP_FILE"

  rm -f ${LOOP_STATUS_FILE}

  stop_check="warning_on_error"

  # check if the tolerance infrastructure needs to be initialized
  if [[ $(hostname) == *"daint"* ]]; then
    echo "Initializing tolerance infrastructure..."
    # make sure that repos do not exist
    if [ -d probtest ]; then
        rm -rf probtest
    fi
    if [ -d icon-test-references ]; then
        rm -rf icon-test-references
    fi
    
    # Clone the testsuite
    git clone git@gitlab.dkrz.de:cscs-sw/probtest.git

    # This variable stores the version of the probtest repository
    PROBTEST_HASH=$(cat scripts/buildbot_scripts/probtest_hash)

    # change to the relevant version of the reference data
    cd probtest
    git reset --hard ${PROBTEST_HASH}
    cd ..

    # Clone the reference data
    git clone git@gitlab.dkrz.de:cscs-sw/icon-test-references.git

    # This variable stores the version of the reference data
    REFERENCE_HASH=$(cat scripts/buildbot_scripts/tolerance_hash)

    # change to the relevant version of the reference data
    cd icon-test-references
    git reset --hard ${REFERENCE_HASH}
    cd ..
  fi

  echo "Run all *.run in run directory"
  cd run
  EXP_FILES=`cat runscript_list`

  case $HOSTNAME in
      rcnl*)
	  batch_system=nqsv
	  slurm_user=`whoami`
       	  typeset -A job_submitted
      ;;
      *)
	  case $submit in
	      sbatch*)
		  batch_system=slurm
		  slurm_user=`whoami`
		  typeset -A job_submitted
	      ;;
	      msub*)
		  batch_system=moab_kit
		  slurm_user=`whoami`
		  typeset -A job_submitted
	      ;;
	      *)
		  batch_system=other
	      ;;
	  esac
      ;;
  esac


  
  for EXP_FILE in $EXP_FILES
  do 
    run_command="$submit ./$EXP_FILE"

    if [ -r $EXP_FILE ]
    then 
        if [[ "other" = "${batch_system}" ]]
        then
          # we assume there is not batch system at all, hence LOG files must be
          # created manually so that buildbot can find them
          run_command="${run_command} > LOG.${EXP_FILE}.o 2>&1"
        fi

        echo "---------------------------------------------------------"
        echo " Submit new Script: ${EXP_FILE} at $(date)"
        echo " $run_command "
        echo "---------------------------------------------------------"

        echo $run_command > submit.$EXP_FILE
        if [ "_$submit" = "_qsub -Wblock=true" ] 
        then
            echo "echo \$? > ${EXP_FILE}.status.2" >> submit.$EXP_FILE
        fi

        chmod +x submit.$EXP_FILE        
        if [[ "$batch_system" = "moab_kit" ]]
        then
          jobid=$(./submit.$EXP_FILE 2>&1 | sed '/^$/d')
        elif [[ "$batch_system" = "nqsv" ]]
        then
          #
          # submit and catch jobID
          echo "jobid=\$(($run_command) 2>&1)" > submit.$EXP_FILE
          echo "jobid=\`echo \$jobid |awk 'match(\$0,/[0-9]+/){print substr(\$0,RSTART,RLENGTH)}'\`" >> submit.$EXP_FILE
          # wait for job to finish and catch error status
          echo "errstat=\$((qwait \${jobid}) 2>&1)" >> submit.$EXP_FILE
          echo "errstat=\`echo \$errstat |awk 'match(\$0,/[0-9]+/){print substr(\$0, RSTART, RLENGTH)}'\`" >> submit.$EXP_FILE
          echo "echo \${errstat} > ${EXP_FILE}.status.2" >> submit.$EXP_FILE

          chmod +x submit.$EXP_FILE
          # submit job
          ./submit.$EXP_FILE &
        else
          ./submit.$EXP_FILE &

    fi


    fi
    if [[ "$batch_system" = "slurm" ]]
    then
        sleep 2
        jobid=`squeue -u ${slurm_user} -h -o '%i' -S '-i' | awk 'NR==1{print $1}'`
        job_submitted["$EXP_FILE"]=$jobid
    fi

    if [[ "$batch_system" = "moab_kit" ]]
    then
        sleep 2
        echo $jobid
        job_submitted["$EXP_FILE"]=$jobid
        job_state="Running"
        while [[ "${job_state}" = "Idle" || "${job_state}" = "Running" ]]
          do
           job_state=`checkjob $jobid  | awk 'NR==4{print $2}'`
          done

    fi
  done

  # wait for all jobs to finish
  wait
  sleep 30

  echo $(pwd)

  
  # print and check the results
  for EXP_FILE in $EXP_FILES
  do
    
    if [ -r $EXP_FILE ]
    then 
      echo "---------------------------------------------------------"
      echo " "
      echo " Start of ${EXP_FILE}"
      echo " "
      echo "---------------------------------------------------------"
      cat LOG.$EXP_FILE.*
      echo "---------------------------------------------------------"
      echo " "
      echo " End of ${EXP_FILE}"
      echo " "
      echo "---------------------------------------------------------"

#       if [[ "$batch_system" = "slurm" ]]
#       then
#           sleep 15
#           slurm_jobid=${job_submitted["$EXP_FILE"]}
#           job_state=COMPLETING
#           while [[ "${job_state}" = "COMPLETING" || "${job_state}" = "RUNNING" ]]
#           do
#               sleep 5
#               job_state=`sacct -j ${slurm_jobid} -n -o jobid,state,exitcode | awk 'NR==1{print $2}'`
#           done
#           case "$job_state" in
#               TIMEOUT)
#                   echo 127 > ${EXP_FILE}.status
#               ;;
#               COMPLETED)
#                   echo 0 > ${EXP_FILE}.status
#               ;;
#               *)
#                   echo 191 > ${EXP_FILE}.status         
#               ;;
#           esac    
#       fi

      exp_status_file=$EXP_FILE.final_status
      if [[ "$batch_system" = "moab_kit" ]]
      then
          sleep 15
          job_state="Running"
          while [[ "${job_state}" = "Idle" || "${job_state}" = "Running" ]]
          do
            job_state=`checkjob $jobid  | awk 'NR==4{print $2}'`
          done
          case "$job_state" in
          Completed*)
                  echo 0 > ${exp_status_file}
              ;;
          esac

      fi 

      STATUS_IN_FILE=255
      if [ -r ${exp_status_file} ]
      then
        STATUS_IN_FILE=`cat ${exp_status_file}`
      else
        if [ -r $EXP_FILE.status.2 ]
        then
          STATUS_IN_FILE=`cat ${EXP_FILE}.status.2`
        fi
      fi

      $stop_check 0 $EXP_FILE
    else
      echo "---------------------------------------------------------"
      echo " ${EXP_FILE} does not exist"
      echo "---------------------------------------------------------"
      echo
    fi
  done
  
  cd ..

  ALL_RUNS_OK=`grep ${TEXT_RUN_FAILED} ${LOOP_STATUS_FILE}`

  if [ $? == 0 ]
  then
    echo "One or more Exp-Runs were not successful"
    RETURN_STATUS=1
  else
    echo "All Exp-Runs were successful"
    RETURN_STATUS=0
  fi



  #==============================================================
  echo "|======================================================|"
  echo "|                                                      |"
  echo "|           Ends                                       |"
  echo "|              $(date)                                 |"
  echo "|                                                      |"
  echo "|======================================================|"
}
#==============================================================================

#------------------------------------------------------------------------------\
# read set-up info
. ./run/set-up.info
submit=$use_sync_submit
#-----------------------------------------------------------------------------
# load ../setting if exists
if [ -a ./setting ]
then
  . ./setting
fi
#-----------------------------------------------------------------------------
run_scripts_submit
#------------------------------------------------------------------------------
# return OK Status
exit ${RETURN_STATUS}

