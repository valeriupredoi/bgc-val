#!/bin/bash
echo "uname -a"
uname -a 
echo "USER: " $USER
#eval "$(ssh-agent)"

echo 'SSH_AUTH_SOCK:' $SSH_AUTH_SOCK  
echo 'SSH_CLIENT:' $SSH_CLIENT 
echo 'SSH_CONNECTION:' $SSH_CONNECTION 
echo 'SSH_TTY:'  $SSH_TTY        
 
python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py


#printf "\n\nssh -A -X mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'\n"
#ssh -X -A mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'

ssh -X -A mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/downloadFromMass.py u-ad980'
ssh -X -A mass-cli1 'ls -lhrt /group_workspaces/jasmin2/ukesm/BGC_data/u-ad980; ln -s  /group_workspaces/jasmin2/ukesm/BGC_data/u-ad980/*1y*grid_[UVWT]* /group_workspaces/jasmin2/ukesm/BGC_data/u-ad980/1y/.'

echo "The end of downloadOnMassCli.sh"

