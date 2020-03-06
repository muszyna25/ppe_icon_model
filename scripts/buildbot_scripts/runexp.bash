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
# Check if the first parameter (return status) is not OK
    echo STATUS_IN_FILE ${STATUS_IN_FILE}
    if [[ $1 -ne 0 || ${STATUS_IN_FILE} -ne 0 ]] 
    then
      if [[ $1 -ne 0 ]] 
      then
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        exit $1
      else
        printf '%-50s : %s\n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        exit ${STATUS_IN_FILE}
      fi
    fi

    printf '%-50s : %s\n' $2 ${TEXT_RUN_OK} >> ../${LOOP_STATUS_FILE}
}

#------------------------------------------------------------------------------------

warning_on_error()
{
# Check if the first parameter (return status) is not OK
    if [[ $1 -ne 0 || ${STATUS_IN_FILE} -ne 0 ]] 
    then
      if [[ $1 -ne 0 ]] 
      then
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        echo "*********** WARNING: script failed  $1 ****************"
      else
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        echo "*********** WARNING: script failed  ${STATUS_IN_FILE} ****************"
      fi
    else
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
# runs the scrpits 
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
  cnt=$(grep -c -i "tolerance" run/runscripts_list)
  if [ $? -gt 1 ]
  then
      echo "Could not find run/runscripts_list ..."
      exit 1
  else
      echo "found tolerance check. Initializing infrastructure..."
  fi
  if [ $cnt -gt 0 ]; then
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
  EXP_FILES=`cat runscripts_list`

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
  
  for EXP_FILE in $EXP_FILES
  do 
    run_command="$submit ./$EXP_FILE"

    if [ -r $EXP_FILE ]
    then 
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

  echo ${pwd}

  
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

