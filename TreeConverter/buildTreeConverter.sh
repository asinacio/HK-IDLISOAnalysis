## S. Jenkins
## Set up everything to build the TreeConverter code
##
## This is a modified version for compiling the SKTreeConverter
## so uses skofl libraries

source /usr/local/sklib_gcc8/skofl_22b/env.sh
source /usr/local/sklib_gcc8/root_v6.22.06_python3.6/bin/thisroot.sh

export ANAUP=../  #/user/sjenkins/HyperK/HK-IDLISOAnalysis

echo "Building TreeConverter."

g++ -Wall SKTreeConverter.cc -o TreeConverter `root-config --cflags --libs` -I $ANAUP/utils

rm ./*~
