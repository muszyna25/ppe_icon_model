#!/bin/ksh

# Define suite base directory
basedir=/home/ms/de/${USER}/ICON_r2B06_30d

# change working directory
cd ${basedir}/sms

# play suite:
cdp << EOF
   myalias
   play icon.def
   sms
   begin icon
EOF



