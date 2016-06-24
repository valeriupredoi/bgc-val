#!/bin/bash
echo "uname -a"
uname -a 
echo "USER: " $USER


echo 'SSH_AUTH_SOCK:' $SSH_AUTH_SOCK  
echo 'SSH_CLIENT:' $SSH_CLIENT 
echo 'SSH_CONNECTION:' $SSH_CONNECTION 
echo 'SSH_TTY:'  $SSH_TTY   

python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py

ssh -X -A jasmin-sci2 'cd /home/users/ldemora/workspace/ukesm-validation; ipython /home/users/ldemora/workspace/ukesm-validation/theWholePackage.py u-ad980 ReportOnly'

rsync -avP /home/users/ldemora/workspace/ukesm-validation/report-u-ad980.tar.gz ledm@pmpc1448.npm.ac.uk:~/ImagesFromJasmin/.
rsync -avP /home/users/ldemora/workspace/ukesm-validation/report-u-ad980.tar.gz ledm@pmpc1486.npm.ac.uk:~/ImagesFromJasmin/.

echo "The end of runReportOnSci1.sh"
