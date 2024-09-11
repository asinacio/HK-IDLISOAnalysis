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
#include "TMinuit.h"

using namespace std;

class Fitters : public TObject{

 public:

  // The constructor and destructor for DataFilter
  Fitters(){ }
  ~Fitters(){ }

  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////

  // Chi2 function
  Double_t Chi2( std::vector<Double_t> wnpe, std::vector<Double_t> d, Double_t intensity, Double_t proposedExL );

  // Log-likelihood function
  Double_t logL();

  // Minimiser: "Manual" Grid Search
  void GridSearch();

  // Minimiser: Minuit
  void Minuit();
  
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
