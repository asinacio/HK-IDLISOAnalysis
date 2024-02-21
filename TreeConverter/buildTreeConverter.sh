## S. Jenkins
## Set up everything to build the TreeConverter code
##
## Some of these paths will be dependent on the structure
## you choose to install WCSim in, and therefore may need
## changing to compile correctly

export LD_LIBRARY_PATH=$WCSIM_HOME/WCSim_build/WCSim_col/:$LD_LIBRARY_PATH
export ANAUP=../  #/user/sjenkins/HyperK/HK-IDLISOAnalysis

echo "Building TreeConverter."

g++ -Wall TreeConverter.cc -o TreeConverter `root-config --cflags --libs` -ltbb -I $WCSIMDIR/include -I $ANAUP/utils -L$WCSIM_HOME/WCSim_build/WCSim_col -lWCSimRoot

rm ./*~
