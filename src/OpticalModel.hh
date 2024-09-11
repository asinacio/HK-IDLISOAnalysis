////////////////////////////////////////////////////////////////////
///
/// FILENAME: OpticalModel.hh
///
/// CLASS: OpticalModel
///
/// BRIEF: A class create the optical model for the fit
///                
/// AUTHOR: Ana Sofia Inacio <ana.carpinteiroinacio@physics.ox.ac.uk>
///
/// DETAIL: This can include different functions for dedicated analyses
///         such as a full optical fit, a scattering fit only, etc
///
////////////////////////////////////////////////////////////////////

#ifndef __OpticalModel__
#define __OpticalModel__

#include "TObject.h"

using namespace std;

class OpticalModel : public TObject{

 public:

  // The constructor and destructor for DataFilter
  OpticalModel(){ }
  ~OpticalModel(){ }

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

  ClassDef(OpticalModel,1)
  
};
#endif
