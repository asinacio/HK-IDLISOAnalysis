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



