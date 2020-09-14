
while [ 1 -eq 1 ];
do
    ./analysis_timeseries.py u-bw717 physics
    ./analysis_timeseries.py u-bx082 physics

    ./analysis_timeseries.py u-bx188 level1;
    #./analysis_timeseries.py u-bw837 level1; 
    ./analysis_timeseries.py u-bv936 level1; 
    ./analysis_timeseries.py  u-bw462 level1;
    #./analysis_timeseries.py u-bv270 physics; 
    #/./analysis_timeseries.py u-bv334 level1; 
#    ./analysis_timeseries.py u-bu847 level1;
#    ./analysis_timeseries.py u-bu794 level1;
    #./analysis_timeseries.py u-bu737 level1; 
#    ./analysis_timeseries.py u-bu504 level1;
    ./analysis_timeseries.py u-bg555 physics;
    ./analysis_timeseries.py u-ar766 physics; 
#/    ./analysis_timeseries.py u-bt931 level1;
#    ./analysis_timeseries.py u-bt670 physics
#/;    ./analysis_timeseries.py u-bt320 level1;
#    ./analysis_timeseries.py u-bt233 level1;
#    ./analysis_timeseries.py u-bt089 physics; 
#/    ./analysis_timeseries.py u-bk093 physics; 
#    ./analysis_timeseries.py u-bs733 physics; 
#    ./analysis_timeseries.py u-bs704 level1; 
#    ./analysis_timeseries.py u-bs522 level1; 
    ./analysis_timeseries.py u-bb446 level1; 
    ./analysis_timeseries.py u-aw310 level1;

    ./analysis_compare.py
    rsync -av /home/users/ldemora/workspace/ukesm-validation/CompareReports /group_workspaces/jasmin4/esmeval/public/.
    echo "sleeping"
    sleep 10000;
done
