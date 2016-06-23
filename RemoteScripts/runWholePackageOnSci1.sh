#!/bin/bash
python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py
ssh -X ldemora@jasmin-sci1 -X 'ipython /home/users/ldemora/workspace/ukesm-validation/theWholePackage.py u-ad980'
