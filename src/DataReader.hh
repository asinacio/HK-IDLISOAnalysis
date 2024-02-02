////////////////////////////////////////////////////////////////////
///
/// FILENAME: DataReader.hh
///
/// CLASS: DataReader
///
/// BRIEF: A class to read input root files containing data for analysis
///                
/// AUTHOR: Ana Sofia Inacio <ana.carpinteiroinacio@physics.ox.ac.uk>
///
/// DETAIL: DataManager reads the root files with structure produced
///         by TreeConverter and stores the data to be used in the fit.
///
////////////////////////////////////////////////////////////////////

#ifndef __DataReader__
#define __DataReader__

#include "TFile.h"
#include "TTree.h"

using namespace std;

class DataReader{

 public:

  // The constructor and destructor for DataManager
  DataReader();
  ~DataReader(){};

  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////

  /////////////////////////////////
  ////////     SETTERS     ////////
  /////////////////////////////////

  /////////////////////////////////
  ////////     GETTERS     ////////
  /////////////////////////////////

 private:


};
#endif
