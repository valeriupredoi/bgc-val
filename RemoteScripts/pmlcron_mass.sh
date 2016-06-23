#!/bin/bash
source /users/modellers/ledm/.bashrc
. ~/.keychain/$HOSTNAME-sh

export XAUTHORITY=/users/modellers/ledm/.Xauthority
export DISPLAY=':0'
#eval `ssh-agent -s` 
#ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh' 
ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/downloadOnMassCli.sh'
#kill $SSH_AGENT_PID

