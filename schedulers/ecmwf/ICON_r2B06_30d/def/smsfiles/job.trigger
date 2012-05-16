#!/bin/ksh
# @ shell    = /bin/ksh
# @ job_name = ICON_TRIGGER
# @ class    = express
##@ class    = ts
# @ job_type = serial
# @ notification = error
# @ notify_user  = Florian.Prill@dwd.de
# @ queue

set -x

date

# trigger ICON SMS suite ------------------------------
cdp << ENDCDP
set SMS_PROG 903409
myalias
force complete /icon/ifs_fct/fct_ifs
status -a
exit
ENDCDP

exit
