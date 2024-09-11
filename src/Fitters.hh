////////////////////////////////////////////////////////////////////
///
/// FILENAME: Fitters.hh
///
/// CLASS: Fitters
///
/// BRIEF: A class which includes fitter methods for the analysis
///                
/// AUTHOR: Ana Sofia Inacio <ana.carpinteiroinacio@physics.ox.ac.uk>
///
/// DETAIL: This can include different fit methods and minimizers to
///         obtain the best fit parameters of the Optical Model
///
////////////////////////////////////////////////////////////////////

#ifndef __Fitters__
#define __Fitters__

#include "TObject.h"

using namespace std;

class Fitters : public TObject{

 public:

  // The constructor and destructor for DataFilter
  Fitters(){ }
  ~Fitters(){ }

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

  ClassDef(Fitters,1)

};
#endif
