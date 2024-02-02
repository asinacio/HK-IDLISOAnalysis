## Set up everything to build the event display code

export LD_LIBRARY_PATH=$WCSIM_HOME/WCSim_col/:$LD_LIBRARY_PATH
export ANAUP=/data/snoplus3/inacio/hyperK/HK-IDLISOAnalysis

echo "Building TreeConverter."

g++ -Wall TreeConverter.cc -o TreeConverter `root-config --cflags --libs` -I $WCSIMDIR/include -I $ANAUP/utils -L$WCSIM_HOME/WCSim_col -L$WCSIM_HOME/WCSimRoot # -ltbb -lWCSimRoot

rm ./*~
