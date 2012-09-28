#!/bin/ksh
# @ shell    = /bin/ksh
# @ job_name = ICON_TRIGGER
# @ class    = express
##@ class    = ts
# @ job_type = serial
# @ notification = error
# @ notify_user  = Florian.Prill@dwd.de
# @ queue

# submitted with
# ecaccess-job-submit -noDirectives -eventIds 167 -queueName ecgate /home/ms/de/dfi0/ICON_r2B06_10d/def/smsfiles/job.trigger 


set -x

date

# trigger ICON SMS suite ------------------------------
cdp << ENDCDP
set SMS_PROG 903409
myalias
force complete /ifs_trigger/fct_ifs
#force queued  /icon/forecast/check_progress
status -a
exit
ENDCDP

exit
