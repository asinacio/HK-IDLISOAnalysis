# Source code for HK-IDLISOAnalysis

Notes for development:

- Logic:
  1) Class DataReader should be responsible for reading the input data files into Run/Source objects and PMT objects
  2) Class DataRunInfo should contain all the run-related/source-related info
  3) (Future) Class DataPMT should contain PMT objects with the corresponding PMT-related info
  4) Class DataFilter should apply PMT/data selection/cuts
  5) Class OpticalModel deals with the optical models for analyses, chi2/likelihood functions, etc
  6) (Future) Class Fit should deal with the minimisation of the model parameters (can have different minimisation methods)
  7) (Future) Some class to "store" the results

*Important*
In order to have any new classes work with the CMake build, a couple of things are required:
  1) Add NewClass.hh to the defined headers list in `src/CMakeLists.txt`
  2) Add NewClass.cc to the `add_library` command in `src/CMakeLists.txt`
  3) Add a line to `LinkDef.h` to link this class to the dictionary built by CMake. Add a new line of `#pragma link C++ class NewClass+`. *Only do this* when the NewClass actually has functions to be added, otherwise `make` won't like it

If there are any issues with this, let me (Sam) know.