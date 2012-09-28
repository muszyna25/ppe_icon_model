#!/bin/ksh

# terminate sms:
cdp << EOF
   myalias
   cancel -fy icon
   halt -y
   terminate -y
EOF

ps -ef | grep sms | grep ${USER}

rm sms.check*
rm icon/*.job*
rm -rf $HOME/sms/*
rm -rf $HOME/sms_server