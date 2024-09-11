#include "Fitters.hh"

ClassImp( Fitters )

//////////////////////////////////////
//////////////////////////////////////

Double_t Fitters::Chi2( std::vector<Double_t> wnpe, std::vector<Double_t> d, Double_t intensity, Double_t proposedExL ){

  Double_t chi2 = 0.0;

  for( Int_t i = 0; i < (Int_t)wnpe.size(); i++ ){
    Double_t model = intensity * exp(-d[i] / proposedExL); // not considering the PMT angular responses yet
    // The model should be loaded from the OpticalModel class!!

    Double_t value = (wnpe[i]-model)*(wnpe[i]-model)/model;

    chi2 = chi2 + value;
  }

  return chi2;

}

//////////////////////////////////////
//////////////////////////////////////

Double_t Fitters::logL(){

  

}

//////////////////////////////////////
//////////////////////////////////////

void Fitters::GridSearch(){

  double minExL = 0.0;
  double maxExL = 10000.0; // in cm (?)
  double bestChi2 = 10e30; // start with absurdly large value for the best Chi2
  double bestExL = 0.0;
  for( int i = minExL; i < maxExL; i++ ){
    double chi2 = CalcChi2(wNPE,d,sourceI,i);
      if(chi2 < bestChi2){
        bestChi2 = chi2;
        bestExL = i;
      }
  }

}

//////////////////////////////////////
//////////////////////////////////////

void Fitters::Minuit(){



}
