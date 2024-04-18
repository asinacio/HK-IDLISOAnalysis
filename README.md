# HK-IDLISOAnalysis
Repository for development of the Hyper-K Laser Injection System Optical Calibration Analysis

## Build Instructions
Requirements:
- C++11 compiler
- CMake 3.5+
- ROOT 6+

Source ROOT before compilation. To install:
```
mkdir build; cd build;
cmake ../
make install
```

To build with support for WCSim, first set up `$WCSIMDIR` to point to the WCSim directory containing `src/` and `lib/`. This is typically set up during WCSim installation. Then `source this_wcsim.sh` from the build. Then:
```
mkdir build; cd	build;
cmake ../ -DUSE_WCSIM=1
make install
```

Finally, return to the top level directory, and source the setup script. This will add the installed executables to your path:
```
cd ../
source setup.sh
```

HK-IDLISOAnalysis is now ready to use!


## TreeConverter

To run:
```

TreeConverter -i [inputfile.root] -o [outputfile.root] -m [macro.mac]
```
Possible arguments:
- `-i` : Input file
- `-o` : Output file
- `-m` : Macro file (the one used to run the input simulation)
- `-b` : Use B&L PMTs only
- `-v` : Print debug info
- `-h` : Display the help message

Only `-i` and `-m` are actually required for the converter to run, everything else will take default values.

